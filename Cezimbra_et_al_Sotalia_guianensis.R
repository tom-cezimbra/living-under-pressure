# Packages ------------------------------------------------------
if(!require(rgdal)) install.packages("rgdal", dependencies = T)
if(!require(sf)) install.packages("sf", dependencies = T)
if(!require(raster)) install.packages("raster", dependencies = T)
if(!require(spThin)) install.packages("spThin", dependencies = T)
if(!require(usdm)) install.packages("usdm", dependencies = T)
if(!require(biomod2)) install.packages("biomod2", dependencies = T)
if(!require(dplyr)) install.packages("dplyr", dependencies = T)
if(!require(ggplot2)) install.packages("ggplot2", dependencies = T)
if(!require(corrplot)) install.packages("corrplot", dependencies = T)
if(!require(parallel)) install.packages("corrplot", dependencies = T)

#### Resetting parallel environment (foreach) ####
env <- foreach:::.foreachGlobals
rm(list = ls(name = env), pos = env)

#### Defining working directory ####
setwd("./Habitat_suitability/Camadas")
Species <- "Sotalia_guianensis"

# Occurrence data
sotalia=read_sf("Merge_sotalia.shp")

setwd("./Habitat_suitability/Camadas/Camadas_krigadas")

#### Environmental variables ####

chl=raster("cloro_new.tif")
depth=raster("Depth_new.tif")
salt=raster("salinidade_new.tif")
slope=raster("Slope_new.tif")
land_dist=raster("coast_new_f.tif")
sotalia_transformado <- st_transform(sotalia, st_crs(chl))

dep = (resample(depth,chl,method="ngb"))
slop = (resample(slope,chl,method="ngb"))
sal = (resample(salt,chl,method="ngb"))


# Merging layers
biostack <- stack(chl, dep, sal, slop)
pal <- colorRampPalette(c("yellow", "blue"))

names(biostack)[1] <- "Chlorophyll"
names(biostack)[2] <- "Depth"
names(biostack)[3] <- "Salinity"
names(biostack)[4] <- "Slope"
#names(biostack)[5] <- "Land_Distance"
plot(biostack, col = pal(20))

setwd("./Habitat_suitability")

dir.create("Output")
dir.create(paste("Output/", Species, sep = ""))
setwd(paste("Output/", Species, sep = ""))

# Correlation tests -------------------------------------------------------

pairs(biostack) # Correlation
biostack_df <- as.data.frame(biostack)
vif(biostack_df) # Multicollinearity
vifstep(biostack_df)
vifcor(biostack_df, th = 0.7)

# Merging layers without collinearity problem
# biostack <- stack(chl, dep, sal, slo)
png("1. Biostack.png", width = 3200, height = 2400, res = 300); par(oma = c(1,1,1,1)); plot(biostack, col = pal(20)); dev.off()

library(dismo)
set.seed(1963) 
backgr <- randomPoints(biostack, 10000)
absclim <- data.frame(raster::extract(biostack, backgr)) # 'extract' function from 'raster' package
absclim.std <- data.frame(scale(absclim)) 
M <- cor(absclim.std,use="pairwise.complete.obs")
corrplot.mixed(M, upper = "ellipse", lower = "number",number.cex = 1.2,tl.cex = 0.5)

# Presence data -----------------------------------------------------------

sot<-cbind(sotalia_transformado$Lat, sotalia_transformado$Long) # Compiling only relevant data (lat, long)
colnames(sot)<-c("Latitude","Longitude")
sot<-as.data.frame(sot) 
sot$Especie="Sotalia guianensis"

sot_na <- sot %>% tidyr::drop_na(Longitude, Latitude)
nrow(sot)-nrow(sot_na)

#Checking for inconsistencies
sot_na$Longitude<-as.numeric(sot_na$Longitude)
sot_na$Latitude<-as.numeric(sot_na$Latitude)
flags_spatial <- CoordinateCleaner::clean_coordinates(
  x = sot_na, 
  species = "Especie",
  lon = "Longitude", 
  lat = "Latitude",
  tests = c("capitals", 
            "centroids", 
            "duplicates", 
            "equal", 
            "gbif", 
            "institutions", 
            "urban", 
            "validity", 
            "zeros" 
  )
)

#Filtered data
sot_f <- sot_na %>% 
  dplyr::filter(flags_spatial$.summary == TRUE)
sot_f = sot_f %>%
  dplyr::select(Especie, Longitude, Latitude)
PresenceData <- sot_f

# Plotting presence data
plot(chl, xlim=c(-44.9,-43.4),ylim=c(-23.4,-22.9), main="Pontos de ocorrencia por Temperatura")
points(PresenceData$Longitude, PresenceData$Latitude, col="red", cex=1)
View(PresenceData)

# Thinning presence data for 2km --------------------------------------------------

dir.create("Thinned")

PresenceDataThinned <- thin( loc.data = PresenceData,
                             lat.col = "Latitude",
                             long.col = "Longitude",
                             spec.col = "Especie",
                             thin.par = 2,
                             reps = 5,
                             locs.thinned.list.return = TRUE,
                             write.files = TRUE,
                             max.files = 5,
                             out.dir = "Thinned",
                             out.base = "Presence_thinned",
                             write.log.file = TRUE,
                             log.file = "Thinned/Presence_thinned_log_file.txt" )

#Loading thinned data
sot_thin_2km=read.csv("Thinned/Presence_thinned_thin1.csv", sep=",")
sot_thin_2km["Occ"] <- c("1")

#### Prepare data & parameters ####

# Select the name of the studied species

myRespName <- "Sotalia_guianensis"

# Get corresponding presence/absence data

myResp <- as.numeric(sot_thin_2km$Occ)
str(myResp)

myRespXY <- sot_thin_2km[which(myResp==1),c("Longitude", "Latitude")]
colnames(myRespXY)


# Pseudo-absences extraction ----------------------------------------------
bm.Species <- BIOMOD_FormatingData(resp.var = myResp,
                                   expl.var = biostack,
                                   resp.xy = myRespXY,
                                   resp.name = myRespName,
                                   PA.nb.rep = 3,
                                   PA.nb.absences = (nrow(sot_thin_2km)*3),
                                   PA.strategy = "disk", 
                                   PA.dist.min = 2000)

# Summary
bm.Species
sum(bm.Species@PA.table[["PA1"]] == TRUE)-nrow(sot_thin_2km) # n of PAs

# Plotting pseudo-absences
plot(bm.Species)
png("2. Pseudo-absences.png", width = 2400, height = 2400, res = 300); plot(bm.Species); dev.off()

#### Parameterize modeling options ####

# Preparing for Maxent

# Saving explanatory variables in .ascii for MaxEnt algorithm

dir.create("maxent_bg")
maxent.background.dat.dir <- "maxent_bg"

for(var_ in names(biostack)){
  cat("/n> saving", paste0(var_, ".asc"))
  writeRaster(subset(biostack, var_), 
              filename = file.path(maxent.background.dat.dir, paste0(var_, ".asc")),
              overwrite = TRUE)
}

# Raster files
caminhos_rasters <- c("./maxent_bg/Chlorophyll.asc",
                      "./maxent_bg/Depth.asc",
                      "./maxent_bg/Salinity.asc",
                      "./maxent_bg/Slope.asc")
                      # Altere e adicione outros caminhos conforme necessário


# Raster objects
rasters <- lapply(caminhos_rasters, raster)

# Check spatial resolution
for (i in seq_along(rasters)) {
  print(paste("Resolução espacial para o raster", i, ":", res(rasters[[i]])))
}

# Check spatial extent
for (i in seq_along(rasters)) {
  print(paste("Extensão espacial para o raster", i, ":", extent(rasters[[i]])))
}

# Check projection
for (i in seq_along(rasters)) {
  print(paste("Projeção para o raster", i, ":", projection(rasters[[i]])))
}

# Correcting incongruences
min_extent <- raster::intersect(rasters[[1]], rasters[[2]], rasters[[3]], rasters[[4]])

library(raster)

rasters <- list(rasters[[1]], rasters[[2]], rasters[[3]], rasters[[4]])
min_extent <- Reduce(raster::intersect, rasters)
rasters_resampled <- lapply(rasters, function(r) {
  projected_raster <- raster::projectRaster(r, min_extent)
  return(projected_raster)
})
for (i in seq_along(rasters_resampled)) {
  print(paste("Extensão espacial para o raster", i, ":", extent(rasters_resampled[[i]])))
}

# Saving resampled layers
for (i in seq_along(rasters_resampled)) {
  writeRaster(rasters_resampled[[i]], 
              filename = file.path(maxent.background.dat.dir, paste0("resampleado_", i, ".asc")), 
              overwrite = TRUE)
}

for (i in seq_along(rasters_resampled)) {
  caminho_salvar <- caminhos_rasters[i]  # Caminho do arquivo original
  
  if (file.exists(caminho_salvar)) {
    file.remove(caminho_salvar)
  }
  
    writeRaster(rasters_resampled[[i]], filename = caminho_salvar, overwrite = TRUE)
}

# Defining the path for maxent.jar file 

path.to.maxent.jar <- file.path("C:/manuscript", "maxent.jar")

# Defining Models Options

bm.opt <- BIOMOD_ModelingOptions(GAM = list(k = 4),
                                 MAXENT = list(path_to_maxent.jar = path.to.maxent.jar,
                                               background_data_dir = maxent.background.dat.dir,
                                               maximumbackground = 10000))

##### Parallel processing setup using doParallel and multicore ####

cl <- makeCluster(2)
registerDoParallel(cl)
on.exit(stopCluster(cl))
options(mc.cores = parallel::detectCores() - 1)

#### Running single models ####

{sotaliamodel <- BIOMOD_Modeling(bm.Species,
                                 modeling.id = paste("model_", Species, sep=""),
                                 models = c("GLM", "GAM", "RF", "GBM", "MAXENT"),
                                 bm.options = bm.opt,
                                 CV.nb.rep = 10,
                                 CV.strategy = "kfold",
                                 CV.k = 5,
                                 CV.perc = 0.7,
                                 do.full.models = F,
                                 prevalence = 0.5,
                                 var.import = 3,
                                 metric.eval = c("TSS", "ROC"))

beepr::beep(8)}

# Get evaluation scores & variables importance
get_variables_importance
ModelEvalsot <- get_evaluations(sotaliamodel)
write.csv(ModelEvalsot, "3. ModelEval_single.csv")

VarImportsot <- get_variables_importance(sotaliamodel)
write.csv(VarImportsot, "4. VarImport_single.csv")

# Represent evaluation scores & variables importance

Graph_EvalSinglesot <- bm_PlotEvalMean(bm.out = sotaliamodel, dataset = "validation")
png("5. EvalSingle.png", width = 2400, height = 2400, res = 300); Graph_EvalSinglesot; dev.off()

Graph_BPEvalSinglesot <- bm_PlotEvalBoxplot(bm.out = sotaliamodel, dataset = "validation", group.by = c('algo', 'algo'))
png("6. BoxplotEvalSingle.png", width = 2400, height = 2400, res = 300); Graph_BPEvalSinglesot; dev.off()

ModelEval_ROC <- ModelEvalsot %>% filter(metric.eval == "ROC")
newalgo <- factor(ModelEval_ROC$algo, levels = c("GLM", "GAM", "RF", "GBM", "MAXENT"))

NewGraph_EvalSinglesot <-  ggplot(ModelEval_ROC, aes(x = newalgo, y = validation, fill = algo)) +
  stat_boxplot(geom = "errorbar") + 
  geom_hline(yintercept = 0.85, linetype='longdash', col='red') +
  labs(title = paste(Species), x = "Algorithms", y = "ROC")+
  guides(fill = "none", color = "none") +
  theme_classic() +
  geom_boxplot() +
  ylim(0.15, 1) +
  theme(legend.title = element_blank(),
        title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_EvalSinglesot)

png("6. NewGraph_EvalSingle.png", width = 1800, height = 1800, res = 300); NewGraph_EvalSinglesot; dev.off()

#### Running ensemble model ####

{sot_EM <- BIOMOD_EnsembleModeling(sotaliamodel,
                                            models.chosen = "all",
                                            em.by = "all",
                                            metric.select = c("ROC"),
                                            em.algo = c("EMmean", "EMwmean", "EMca", "EMcv"),
                                            metric.select.thresh = c(0.85),
                                            metric.eval = c("ROC"),
                                            var.import = 3)

beepr::beep(8)}

summary(sot_EM)
get_predictions(sot_EM)

{sot_EM1 <- BIOMOD_EnsembleModeling(sotaliamodel,
                                   models.chosen = "all",
                                   em.by = "PA+run",
                                   metric.select = c("ROC"),
                                   em.algo = c("EMwmean"),
                                   metric.select.thresh = c(0.85),
                                   metric.eval = c("ROC"),
                                   var.import = 3)
  
  beepr::beep(8)}

summary(sot_EM1)
get_predictions(sot_EM1)
#### Evaluations for the ensemble model ####
evalsotensemble1 = get_evaluations(sot_EM1)

#### Plotting evaluations for the ensemble model ####
ens_sot_eval =  bm_PlotEvalMean(bm.out = sot_EM)

length(sot_EM)
summary(sot_EM@em.models_kept)

# Get evaluation scores & variables importance

ModelEval_EM <- get_evaluations(sot_EM)
write.csv(ModelEval_EM,paste0("7. ModelEval_EM.csv"))

ModelEval_EM1 <- get_evaluations(sot_EM1)
write.csv(ModelEval_EM1,paste0("7.1. ModelEval_EM1.csv"))

VarImport_EM <- get_variables_importance(sot_EM)
write.csv(VarImport_EM, "8. VarImport_EM.csv")

# Represent evaluation scores & variables importance

Graph_BoxplotVarAlgo <- bm_PlotVarImpBoxplot(bm.out = sot_EM, group.by = c('expl.var', 'algo', 'algo'))
png("9. BoxplotVarAlgo.png", width = 2400, height = 2400, res = 300); Graph_BoxplotVarAlgo; dev.off()
?bm_PlotVarImpBoxplot()

Graph_BoxplotAlgoVar <- bm_PlotVarImpBoxplot(bm.out = sot_EM, group.by = c('algo', 'expl.var', 'algo'))
png("10. BoxplotAlgoVar.png", width = 2400, height = 2400, res = 300); Graph_BoxplotAlgoVar; dev.off()

NewGraph_BoxplotVarImp <-  Graph_BoxplotVarAlgo$tab %>%
  filter(algo == "EMwmean") %>%
  ggplot() +
  geom_boxplot(aes(x = expl.var, y = var.imp, fill = factor(expl.var), color = factor(expl.var))) +
  scale_fill_manual(values = c("white", "white", "white", "white", "white", "white")) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "cyan3", "brown4")) +
  theme_classic() +
  labs(title = paste(Species), x = "Explanatory variables", y = "Variable importance") +
  guides(fill = "none", color = "none") +
  ylim(0, 1) +
  theme(title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_BoxplotVarImp)

png("11. NewBoxplot_wmean.png", width = 1800, height = 1800, res = 300); NewGraph_BoxplotVarImp; dev.off()

# Represent response curves
Graph_ResponseCurves_wmean <- bm_PlotResponseCurves(bm.out = sot_EM,
                                                    models.chosen = get_built_models(sot_EM)[2],
                                                    fixed.var = 'median')
png("12. ResponseCurves_wmean.png", width = 2400, height = 2400, res = 300); Graph_ResponseCurves_wmean; dev.off()

ggdat <- Graph_ResponseCurves_wmean$tab
ggdat$algo <- ifelse(grepl("_GLM", ggdat$pred.name, ignore.case = T), "GLM",
                     ifelse(grepl("_GAM", ggdat$pred.name, ignore.case = T ), "GAM",
                            ifelse(grepl("_GBM", ggdat$pred.name, ignore.case = T ), "GBM",
                                   ifelse(grepl("_RF", ggdat$pred.name, ignore.case = T ), "RF",
                                          ifelse(grepl("_MAXENT", ggdat$pred.name, ignore.case = T ), "MAXENT", "NA")))))

ggdat$algo <- as.factor(ggdat$algo)

ggdat$id <- as.factor(ggdat$id)

NewResponseCurve <- ggplot(ggdat, aes(x = expl.val, y = pred.val)) +
  geom_smooth(aes(group = algo, color = algo), fill = "indianred1", se = TRUE, method = "gam", alpha = 0.01) +
  geom_smooth(color = "red", fill = "indianred1", se = TRUE, method = "gam") +
  facet_wrap(~ expl.name, scales = "free_x") +
  labs(title = paste(Species), x = "Explanatory variables value", y = "Suitability value") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "",
        title = element_text(size = 16, color = "gray33"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "gray33"),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(1, "lines"))

plot(NewResponseCurve)

png("13. NewResponseCurves_wmean.png", width = 2400, height = 1800, res = 300); NewResponseCurve; dev.off()

#### Projecting single models ####

{Proj_sot <- BIOMOD_Projection(bm.mod = sotaliamodel,
                                        proj.name = "Current",
                                        new.env = biostack,
                                        models.chosen = "all",
                                        metric.binary = "all",
                                        metric.filter = "all",
                                        build.clamping.mask = F,
                                        keep.in.memory = F)

beepr::beep(8)}


#### Projecting ensemble model ####

{Proj_sot_EM <- BIOMOD_EnsembleForecasting(bm.em = sot_EM,
                                                    bm.proj = Proj_sot,
                                                    proj.name = "Current_EM",
                                                    models.chosen = "all",
                                                    metric.binary = "all",
                                                    metric.filter = "all",
                                                    output.format = ".img")

beepr::beep(8)}

Proj_sot_EM@modeling.id

save.image(file = (paste("../../", Species, ".RData", sep = "")))

Model_Eval2km<- read.csv(file = "C:/manuscript/7. ModelEval_EM.csv") 
pesos<-read.csv(file = "C:/manuscript/final models weights.csv",col.names = "pesos")
frequencia <- table(pesos)
print(frequencia)


# PA+RUN ---------------------------------------------------


#### Evaluations for the ensemble model ####
evalsotensemble1 = get_evaluations(sot_EM1)
# Get evaluation scores & variables importance

ModelEval_EM1 <- get_evaluations(sot_EM1)
write.csv(ModelEval_EM1,paste0("7.1. ModelEval_EM1.csv"))

VarImport_EM1 <- get_variables_importance(sot_EM1)
write.csv(VarImport_EM1, "8.1. VarImport_EM1.csv")

# Represent evaluation scores & variables importance

Graph_BoxplotVarAlgo1 <- bm_PlotVarImpBoxplot(bm.out = sot_EM1, group.by = c('expl.var', 'algo', 'algo'))
png("9.1. BoxplotVarAlgo1.png", width = 2400, height = 2400, res = 300); Graph_BoxplotVarAlgo1; dev.off()
?bm_PlotVarImpBoxplot()

Graph_BoxplotAlgoVar1 <- bm_PlotVarImpBoxplot(bm.out = sot_EM1, group.by = c('algo', 'expl.var', 'algo'))
png("10.1. BoxplotAlgoVar1.png", width = 2400, height = 2400, res = 300); Graph_BoxplotAlgoVar1; dev.off()

NewGraph_BoxplotVarImp1 <-  Graph_BoxplotVarAlgo1$tab %>%
  filter(algo == "EMwmean") %>%
  ggplot() +
  geom_boxplot(aes(x = expl.var, y = var.imp, fill = factor(expl.var), color = factor(expl.var))) +
  scale_fill_manual(values = c("white", "white", "white", "white", "white", "white")) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "cyan3", "brown4")) +
  theme_classic() +
  labs(title = paste(Species), x = "Explanatory variables", y = "Variable importance") +
  guides(fill = "none", color = "none") +
  ylim(0, 1) +
  theme(title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_BoxplotVarImp1)

png("11.1. NewBoxplot_wmean1.png", width = 1800, height = 1800, res = 300); NewGraph_BoxplotVarImp1; dev.off()

# Represent response curves

Graph_ResponseCurves_wmean1 <- bm_PlotResponseCurves(bm.out = sot_EM1,
                                                    models.chosen = get_built_models(sot_EM1)[1],
                                                    fixed.var = 'median')
png("12.1. ResponseCurves_wmean1b.png", width = 2400, height = 2400, res = 300); Graph_ResponseCurves_wmean1; dev.off()

ggdat1 <- Graph_ResponseCurves_wmean1$tab
ggdat1$algo <- ifelse(grepl("_GLM", ggdat1$pred.name, ignore.case = T), "GLM",
                     ifelse(grepl("_GAM", ggdat1$pred.name, ignore.case = T ), "GAM",
                            ifelse(grepl("_GBM", ggdat1$pred.name, ignore.case = T ), "GBM",
                                   ifelse(grepl("_RF", ggdat1$pred.name, ignore.case = T ), "RF",
                                          ifelse(grepl("_MAXENT", ggdat1$pred.name, ignore.case = T ), "MAXENT", "NA")))))

ggdat1$algo <- as.factor(ggdat1$algo)

ggdat1$id <- as.factor(ggdat1$id)

NewResponseCurve1 <- ggplot(ggdat1, aes(x = expl.val, y = pred.val)) +
  geom_smooth(aes(group = algo, color = algo), fill = "indianred1", se = TRUE, method = "gam", alpha = 0.01) +
  geom_smooth(color = "red", fill = "indianred1", se = TRUE, method = "gam") +
  facet_wrap(~ expl.name, scales = "free_x") +
  labs(title = paste(Species), x = "Explanatory variables value", y = "Suitability value") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "",
        title = element_text(size = 16, color = "gray33"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "gray33"),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(1, "lines"))

plot(NewResponseCurve1)

png("13.1. NewResponseCurves_wmean1.png", width = 2400, height = 1800, res = 300); NewResponseCurve1; dev.off()

#### Projecting single models ####

{Proj_sot <- BIOMOD_Projection(bm.mod = sotaliamodel,
                               proj.name = "Current",
                               new.env = biostack,
                               models.chosen = "all",
                               metric.binary = "all",
                               metric.filter = "all",
                               build.clamping.mask = F,
                               keep.in.memory = F)

beepr::beep(8)}


#### Projecting ensemble model ####

{Proj_sot_EM <- BIOMOD_EnsembleForecasting(bm.em = sot_EM,
                                           bm.proj = Proj_sot,
                                           proj.name = "Current_EM",
                                           models.chosen = "all",
                                           metric.binary = "all",
                                           metric.filter = "all",
                                           output.format = ".img")

beepr::beep(8)}

Proj_sot_EM@modeling.id

save.image(file = (paste("../../", Species, ".RData", sep = "")))
#### The end :) ####