
# For everything else
library(maxnet)
library(randomForest)
library(gam)
library(terra)
library(PresenceAbsence)
library(dismo)
library(usdm)
library(mapview)
library(sf)
library(biomod2)
library(car)
library(lattice)
library(ggtext)
library(rasterVis)
library(plyr)
library(gridExtra)#for displaying ggplots in one window
library(dplyr)
library(sp)
library(ecospat)
library(biomod2)
library(tidyr)
library(tidyterra)
library(tidyverse)




#  Action points in this coding chunk: 

#setwd("/Users/catherinemercer/Documents/Taiwan_SDMs/SDM_Feb18/Presence_Data") 

inventory <- read_csv("/Users/catherinemercer/Documents/Taiwan_SDMs/SDM_Feb18/Presence_Data/NFI_Inventory_Named_Cleaned.csv")%>%
  dplyr::select(-...1)

names(inventory)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # #  Formating inventory: Adult Trees and Juvenile Trees # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# Cleaning also found in the r script:  

nfi_data <- inventory
# Data cleaning: Remove outliers, and filter based on given criteria
cleaned_data <- nfi_data %>%
  filter(!is.na(Elev1) & !is.na(DBH) & !is.na(Plot_X) & !is.na(Plot_Y)) %>%
  filter(DBH <= quantile(DBH, 0.9999), Plot_X <= 8.8, Plot_Y <= quantile(Plot_Y, 0.9999)) %>%
  group_by(Plot_no) %>%
  filter(sum(BasalArea_m2) >= 0.2 & sum(BasalArea_m2) <= quantile(BasalArea_m2, 0.9999)) %>%
  ungroup()
# Removing plantation plots based on basal area > 90% for target species
plantation_species <- c("Taiwania cryptomerioides", "Cryptomeria japonica", "Chamaecyparis formosensis")
cleaned_data <- nfi_data %>%
  group_by(Plot_no) %>%
  mutate(basal_area_percentage = sum(BasalArea_m2[Sci_Name %in% plantation_species]) / sum(BasalArea_m2)) %>%
  ungroup() %>%
  filter(basal_area_percentage < 0.9)

 
nfi_data <- cleaned_data


# Define life stages: Juveniles and Adults based on DBH percentiles
life_stage_1<- nfi_data %>%
  group_by(Sci_Name) %>%
  mutate(life_stage = case_when(
    DBH < quantile(DBH, 0.25) ~ "juvenile",
    DBH > quantile(DBH, 0.50) ~ "adult",
    TRUE ~ "intermediate"
  )) %>%
  filter(life_stage %in% c("juvenile", "adult"))

table(life_stage_1$Sci_Name) %>% #  how many plots the species is found in 
  sort()

# Seperate inventory data so that you're only moving forwards with species over 50 presences
adults_in_50_plots <- life_stage_1 %>%
  filter(life_stage == "adult") %>%       # keep only adults
  group_by(Sci_Name) %>%                  # group by species
  summarise(n_plots = n_distinct(Plot_no)) %>%  # count unique plots
  filter(n_plots >= 50) %>%              # keep species in ≥50 plots
  pull(Sci_Name) 

# Make sure using the same species for adult and juvenile data  
species_in_50_plots <-  life_stage_1 %>%
  filter(Sci_Name %in% adults_in_50_plots)

# ADD SRTM ELEVATION 

print(species_in_50_plots) # Adults and juveniles

file_path_Elev <- "/Users/catherinemercer/Downloads/elevation.tif"
Elev <- rast(file_path_Elev)
species_in_50_plots$elev_srtm <-terra::extract(Elev,species_in_50_plots[,c('Plot_X_dec','Plot_Y_dec')], ID=FALSE)[,1] 


# SEPERATE ADULTS FOR THE BIOMOD2 PROCESSING 
species_in_50_plots

adult_inventory <-  species_in_50_plots %>%
  filter(life_stage=="adult")

# SUMMARISE BY PLOT
SDM_inventory <- adult_inventory %>%
                  group_by(Plot_no, life_stage)  %>%
                  dplyr::summarise(num_trees = n(), Long=mean(Plot_X_dec), Lat=mean(Plot_Y_dec), Elevation=mean(elev_srtm))

# ABUNDANCE -  Calculate number of trees for each species per plot

species_abund_df <- as.data.frame(table(adult_inventory$Plot_no, adult_inventory$Sci_Name))  # Convert matrix to data frame
species_abund_reshaped <- spread(species_abund_df, key = Var2, value = Freq, fill = 0)
colnames(species_abund_reshaped)[1] <- "Plot_no"
SDM_inventory_abund <- merge(SDM_inventory, species_abund_reshaped, by.x = "Plot_no", by.y = "Plot_no", all.x = TRUE)

SDM_inventory_abund

filepath  <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Post_SDM_Analysis/SDM_inventory_abund.csv"
write_csv(SDM_inventory_abund, file=filepath)
#Export dataframe
# 🌸🌸🌸🌸🌸🌸  move to cleaned_analysis_130326.R 🌸🌸🌸🌸🌸🌸


species_names <- as.list(adult_inventory$Sci_Name)

for (species in species_names) {#cat(species, "\n")   # Print species name
  # Extract presence/absence data for the species
  presence_absence <- ifelse(SDM_inventory_abund[, species] > 0, 1, 0)  
  # Assign the presence/absence data to the corresponding column in SDM_inventory
  SDM_inventory[[species]] <- presence_absence
  # Display the table of 0's and 1's in the column
  #  print(table(SDM_inventory[[species]]))
  # Print an empty line for better readability
  # cat("\n")
}

summary(SDM_inventory$sp)
nrow(SDM_inventory)

SDM_inventory_abv50<- SDM_inventory
nrow(SDM_inventory_abv50)
#filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/111025/SDM_inventory_abv50_100226.csv"
#write_csv(SDM_inventory_abv50, file=filename)


#//////////////////////////////////////////////////////////////////////////////// 
#////////////////// START  ////////////////////////// 
#//////////////////////////////////////////////////////////////////////////////// 

SDM_inventory_abv50 <- read.csv("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/111025/SDM_inventory_abv50_100226.csv", row.names=1)

SDM_inventory_abv50$Plot_no<- row.names(SDM_inventory_abv50)
nrow(SDM_inventory_abv50)
#//////////////////////////////////////////////////////////////////////////////// 
#////////////////// Current: Climate and topography ////////////////////////// 
#//////////////////////////////////////////////////////////////////////////////// 

#Sorting the climate APdata (long to wide format + Left join)

# SUMMARISING CURRENTCLIMATE DATA FROM CLIMATEAP

# currentclimate_1981_1990 <-read.csv("/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/Input_Decade_1981_1990Y.csv")
# currentclimate_1991_2000 <-read.csv("/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/Input_Decade_1991_2000Y.csv")
# currentclimate_2001_2010 <-read.csv("/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/Input_Decade_2001_2010Y.csv")

# names(currentclimate_1981_1990)

# ensamb <- bind_rows(currentclimate_1981_1990, currentclimate_1991_2000, currentclimate_2001_2010)

# filename <- "/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/CombinedCurrentClimate.csv"
# write_csv(ensamb, file=filename)

ensamb <- read.csv("/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/CombinedCurrentClimate.csv")


# SUMMARISE BY PLOT
Current_climate <- ensamb %>%
  group_by(ID1)  %>%
  dplyr::summarise(Long=mean(Longitude), Lat=mean(Latitude), Elevation=mean(Elevation), MAT=mean(MAT), MWMT=mean(MWMT), MCMT=mean(MCMT), TD=mean(TD), MAP=mean(MAP),AHM=mean(AHM),DD_0=mean(DD_0),DD5=mean(DD5), DD_18=mean(DD_18), DD18=mean(DD18),NFFD=mean(NFFD), PAS=mean(PAS), EMT=mean(EMT), EXT=mean(EXT),Eref=mean(Eref),  CMD=mean(CMD), RH=mean(RH))


SDM_inventory_abv50 <- left_join(SDM_inventory_abv50,Current_climate)
SDM_inventory_abv50 <- na.omit(SDM_inventory_abv50)



#//////////////////////////////////////////////////////////////////////////////// 
#////////////////// Topography ////////////////////////// 
#//////////////////////////////////////////////////////////////////////////////// 



file_path_Slope <- "/Users/catherinemercer/Downloads/slope.tif"
Slope <- rast(file_path_Slope)
Slope_csv<-terra::extract(Slope,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE) 
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, Slope_csv)

# file_path_Elev <- "/Users/catherinemercer/Downloads/elevation.tif"
# Elev <- rast(file_path_Elev)
# Elev_csv<-terra::extract(Elev,SDM_inventory_abv50[,c('Lat','Long')], ID=FALSE) 
# SDM_inventory_abv50<-cbind(SDM_inventory_abv50, Elev_csv)
# Aspect 
file_path_Aspect<- "/Users/catherinemercer/Downloads/aspect.tif"
Aspect <- rast(file_path_Aspect)

aspect_radians <- Aspect * pi / 180
# Eastness (sin of aspect)
eastness <- sin(aspect_radians)
Eastness_csv<-terra::extract(eastness,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE)  %>%  
  rename(eastness = "aspect")
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, Eastness_csv)

# Northness (cos of aspect)
northness <- cos(aspect_radians)
northness_csv<-terra::extract(northness,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE)  %>%  
  rename(northness = "aspect")
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, northness_csv)

TPI <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/Topographic/TopographicPositionIndex.tif"
TPI <- rast(TPI)
TPI_csv<-terra::extract(TPI,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE) 
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, TPI_csv)


TWI_clean <- rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/TWI_renameNA_inR_181125.tif")
plot(TWI_clean)
TWI_csv <- terra::extract(TWI_clean, SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE)
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, TWI_csv)

#plot(TPI)
crs(northness)




curv_vert <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/Topographic/Curvature/VerticalCurvature_twnclip_280925.tif"
curv_vert <- rast(curv_vert)
curv_vert_csv<-terra::extract(curv_vert,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE) 
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, curv_vert_csv)

curv_horiz <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/Topographic/Curvature/HorizontalCurvature_twnclip_280925.tif"
curv_horiz <- rast(curv_horiz)
curv_horiz_csv<-terra::extract(curv_horiz,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE) 
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, curv_horiz_csv)

curve_mean <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/Topographic/Curvature/MeanCurvature_twnclip_280925.tif"
curve_mean <- rast(curve_mean)
curve_mean_csv<-terra::extract(curve_mean,SDM_inventory_abv50[,c('Long','Lat')], ID=FALSE) 
SDM_inventory_abv50<-cbind(SDM_inventory_abv50, curve_mean_csv)



#//////////////////////////////////////////////// 
#///////////////// VIF //////////////// 
#//////////////////////////////////////////////// 

    #Clim <- SDM_inventory_abv50 %>%
    #select(MWMT, MCMT, TD, MAP, AHM, DD5, CMD, slope:MeanCurvature)
      # usdm::vif(Clim)
      #summary(Clim)
    
#SDM_inventory_abv50 <- na.omit(SDM_inventory_abv50)

# 🌸🌸🌸🌸🌸🌸🌸🌸🌸🌸 Export here after ommitting NA's 🌸🌸🌸🌸🌸🌸🌸🌸🌸🌸

#nrow(SDM_inventory_abv50)
#row.names(SDM_inventory_abv50) <-SDM_inventory_abv50$Plot_no
#filepath <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Post_SDM_Analysis/SDM_inventory_naomit_110326.csv"
#write_csv(SDM_inventory_abv50, file=filepath)
#nrow(SDM_inventory_abv50)

# 🌸🌸🌸🌸🌸🌸🌸🌸🌸🌸 Continue with VIF 🌸🌸🌸🌸🌸🌸🌸🌸🌸🌸

    #vif_output <- vifstep(Clim)
    #cor(Clim)
    #vif_output <- vifstep(Clim)
    
    # nrow(SDM_inventory_abv50)
    
    Clim <- SDM_inventory_abv50 %>%
      select(DD5, TD, MAP, CMD, slope, eastness, northness, TopographicPositionIndex, TWI_4326_181125.tif,MeanCurvature) %>% rename(TWI = "TWI_4326_181125.tif")
    row.names(Clim) <- SDM_inventory_abv50$Plot_no
    

    # m <- lm(MWMT ~ MCMT + TD + MAP + AHM + DD5 +CMD+slope+eastness+northness+TopographicPositionIndex+TWI_4326_181125.tif+VerticalCurvature+HorizontalCurvature+MeanCurvature, data = Clim)
    
    #library(performance)
    #x <- check_collinearity(m)
    #> # Check for Multicollinearit
    #> x <- check_collinearity(m)
    
    
    
    #barplot(Clim, col = "skyblue", main = "Variance Inflation Factor (VIF)")
    
  

#/////////////////////////////////////////////////////////////////////////////
#////////////////////////////// BIOMOD2 WORKFLOW ////////////////////////
    
 
# Species names for loop WILL NEED CHANGED 
sp_names= names(SDM_inventory_abv50)[(6:86)]

# Remove this file so new one added below 
#file.remove("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/231025/Ensemble_forecast_out/fut_projection_df.csv")
#which(sp_names == "Zelkova.serrata")

#for(sp in sp_names[-(1:which(sp_names == "Zelkova.serrata"))]){ # If loop stops 
setwd("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Models")

print('Have you deleted the files if re-running?')
#start_species <- match("Cinnamomum.subavenium", sp_names)
#for (sp in sp_names[start_species:length(sp_names)]) { 
for (sp in sp_names[64:81]) { # from"Diospyros.morrisiana"
  # Deutzia.pulchra  Daphniphyllum.glaucescens
  # Daphniphyllum.glaucescens not working, run again> 
  # Daphniphyllum.glaucescens not working, remove for now
  # from"Diospyros.morrisiana"
#for (sp in sp_names) { 
print(sp)
  sp <- "Abies.kawakamii"
  
  #//////////////  Presence-absence Extract  ///////////////////# 
  
  PA <- SDM_inventory_abv50  %>% dplyr::select(sp,Long, Lat)
  Vars <- Clim
  #class(Vars)
  #Vars<-terra::extract(Clim, PA[,c('x','y')], ID=FALSE) # extracting the climate data here to an excel file - replace for individual points, addingn them to a dataframe (the rasters just made it easier to extract )
  SpData<-cbind(PA, Vars)
  SpData<-SpData[complete.cases(SpData),] 
  nrow(SpData)
  SpData[]
  library("dismo")
  allModels <- c("GLM","GAM",  "ANN", "MAXENT")  #RF

  dataset <- BIOMOD_FormatingData(
    resp.var = SpData[[sp]], #response variable (0, 1)
    expl.var = Clim, #explanatory variables provided as a raster stack
    resp.name = sp, #name of the species
    resp.xy = SpData[,c('Long','Lat')], #coordinates of the response variable
  )
  SpData[]
  
#//////////////////////////////////////////////////////
#////////////// Tuning the models  ///////////////////
#////////////////////////////////////////////////////
  library(biomod2)
  opt.d <- bm_ModelingOptions(data.type = 'binary',
                              models = allModels,
                              strategy = 'default',
                              bm.format=dataset)
  print("Gam Warning okay")
  opt.d@options$MAXENT.binary.MAXENT.MAXENT@args.values$`_allData_allRun`$path_to_maxent.jar="/Users/catherinemercer/Downloads/maxent"
  
  # This should be doing three repeats of 80/20 splits for the tuning below
  ctrl.train <- caret::trainControl(method = "LGOCV", number=3, p=0.8,
                                    summaryFunction = caret::twoClassSummary, classProbs = TRUE, returnData = FALSE, allowParallel = FALSE)
  

  
  library(doParallel)
  
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  getDoParWorkers()
  
  tuned.ann <- bm_Tuning(model = 'ANN',
                         tuning.fun = 'avNNet', ## see in ModelsTable
                         do.formula = FALSE,
                         bm.options = opt.d@options$ANN.binary.nnet.nnet,
                         bm.format = dataset,
                         ctrl.train = ctrl.train)
  
  tuned.maxent <- bm_Tuning(model = 'MAXENT', ## see in ModelsTable
                            tuning.fun = 'ENMevaluate', ## see in ModelsTable
                            do.formula = TRUE,
                            bm.options = opt.d@options$MAXENT.binary.MAXENT.MAXENT,
                            bm.format = dataset,
                            ctrl.train = ctrl.train,
                            metric.eval = "AICc" ) ##default is TSS)
  
  

  opt.d@options$ANN.binary.nnet.nnet@args.values$`_allData_allRun` <- tuned.ann$`_allData_allRun`
  opt.d@options$MAXENT.binary.MAXENT.MAXENT@args.values$`_allData_allRun`<- tuned.maxent$`_allData_allRun`
  opt.d@options$MAXENT.binary.MAXENT.MAXENT@args.values$`_allData_allRun`$path_to_maxent.jar="/Users/catherinemercer/Downloads/maxent"

  

  unregister_dopar <- function() {  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)}
  unregister_dopar()
  
  print("Tuning Done")

  
  
#///////////////////////////////////////////////////
#////////////// BUILD MODELS - CURRENT CLIMATE    /////
#///////////////////////////////////////////////////
  
citation()
citation("biomod2")
  N_Runs=10  # 10 to match the number of environmental variables (Hill et al 2016)
  
  models <- BIOMOD_Modeling(#GAM = list(type = 's_smoother'), 
                              bm.format=dataset,
                              models = c("GLM","GAM", "ANN", "MAXENT"),  # add in others to justify choice
                              # bm.options = opt.d, 
                              OPT.user = opt.d,
                              OPT.strategy = "user.defined",
                              CV.nb.rep = 5, # increase to 5 - smallest value to redo ,
                              CV.perc = 0.8, #proportion of data used in the training dataset
                              var.import = N_Runs, # default 3 - number of runs to estimate variable importance
                              metric.eval = c("TSS","ROC", "KAPPA" ), 
                              CV.do.full.models = FALSE,
                              scale.models = TRUE,
                              modeling.id = paste(sp,"FirstModeling",sep=""))
    print("Warnings okay. Including 'This error ok: Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred'")
  
  # Plots: evaluation and response curves

  folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Bmd_model_EvalPlots"
              pdf_file <- file.path(folder_path, paste0("bm_model_evaluations_", sp, ".pdf"))
              pdf(pdf_file)
              bm_PlotEvalMean(models)
              bm_PlotEvalBoxplot(models, group.by = c('algo', 'PA'), dataset="validation")# better plots using the data=validation 
              bm_PlotVarImpBoxplot(models,  group.by = c('expl.var','algo' ,'PA'))
              RespCurve<- bm_PlotResponseCurves(models)
              RespCurve$tab# saving response curves
              dev.off() # Close the pdf file
  
              

  # Export response curves data for Ggplot plotting
  filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Bmd_model_EvalPlots/responsecurves.csv"
              write_csv(RespCurve$tab, file=filename, append = sp!=sp_names[1])
  
 
  #Variable importance 
  imp<-get_variables_importance(models)
              imp
              imp_df <- as.data.frame(imp)
              imp_df$sp=sp
              folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Bmd_model_VarImport"
              filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Bmd_model_VarImport/Variable_importance.csv"
              write_csv(imp_df, file=filename, append = sp!=sp_names[1])
  
 
  eval_models <- get_evaluations(models)
              print(eval_models)
              eval_models_df <- as.data.frame(eval_models)
              eval_models_df$sp=sp
              folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Bmd_model_EvalMetrics"
              filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Bmd_model_EvalMetrics/Eval_metrics.csv"
              write_csv(eval_models_df, file=filename, append = sp!=sp_names[1])
  
print("Got past model building")
              
#///////////////////////////////////////////////////
#//////////////  BUILD ENSEMBLE MODEL PREDICTIONS  /////
#///////////////////////////////////////////////////
  
  
  myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = models,
                                        models.chosen = 'all',
                                        em.by = 'PA+run',
                                        em.algo = c('EMwmean'),
                                        metric.select = c('TSS'),
                                        metric.select.thresh = c(0),
                                        metric.select.dataset = 'validation',
                                        metric.eval = c('TSS','ROC','KAPPA'),
                                        var.import = 3,
                                        EMci.alpha = 0.05,
                                        EMwmean.decay = 'proportional')
  
  print("Got past ensemble building")
  
  # EVALUATION AND RESPONE CURVE PLOTS
  folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_EvalPlots"
              pdf_file <- file.path(folder_path, paste0("Ensemble_Bmd_model_EvalPlots", sp, ".pdf"))
              pdf(pdf_file)
              #bm_PlotEvalMean(myBiomodEM)
              # group.by[1] must be 'full.name', 'merged.by.PA', 'merged.by.run' or 'algo'
              bm_PlotEvalBoxplot(myBiomodEM, group.by = c('algo', 'merged.by.PA'))# better plots using the data=validation, tuning of this function different from ensemble
              bm_PlotVarImpBoxplot(myBiomodEM,  group.by = c('algo' ,'merged.by.PA', 'expl.var'))
              RespCurve<- bm_PlotResponseCurves(myBiomodEM)
              RespCurve$tab# saving response curves
              dev.off() # Close the pdf file
            
  # Export 
  filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_EvalPlots/responsecurves.csv"
              write_csv(RespCurve$tab, file=filename, append = sp!=sp_names[1])

  #Variable importance 
  imp<-get_variables_importance(myBiomodEM)
              imp_df <- as.data.frame(imp)
              imp_df$sp=sp
              folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_VarImport"
              filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_VarImport/EnsembleVariable_importance.csv"
              write_csv(imp_df, file=filename, append = sp!=sp_names[1])


  eval_models <- get_evaluations(myBiomodEM)
              eval_models_df <- as.data.frame(eval_models)
              eval_models_df$sp=sp
              folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_EvalMetrics"
              filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_EvalMetrics/Ensemble_Eval_metrics.csv"
              write_csv(eval_models_df, file=filename, append = sp!=sp_names[1])
              
              eval_models_summary <- eval_models %>%
                                      group_by(metric.eval) %>%
                                      summarise(cutoff=mean(cutoff), sensitivity = mean(sensitivity), specificity= mean(specificity), calibration = mean(calibration), validation = mean(validation)) %>%
                                      as.data.frame() %>%
                                      mutate(sp=sp)
              folder_path <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_EvalMetrics"
              filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_EvalMetrics/Ensemble_eval_models_summary.csv"
              write_csv(eval_models_summary, file=filename, append = sp!=sp_names[1])
  
              
              
              
#///////////////////////////////////////////////////
#////////////// PREDICT ONTO CURRENT CLIMATE (individual models and ensemble)  /////
#///////////////////////////////////////////////////          

  print('made it to projecting current climate')
  #//////////////  Single models  ///////////////////# 
  
  preds_pres <- BIOMOD_Projection(  bm.mod = models,
                                      new.env = Clim,
                                      proj.name = 'CurrentRangeSize',
                                      selected.models = 'all',
                                      binary.meth = 'TSS' # Finds the optimal TSS threshold that maximises the tss threshold, applies the threshols and saves binary raster 
   )
  
  # preds_pres@proj.out@val  
  # get

  # binary_maps <- get_predictions(preds_pres, as.data.frame = FALSE, binary = TRUE)
              # binary = TRUE ensures you extract the binary version (presence/absence) of the maps.
              # The result will be a RasterStack (or SpatRaster if you use terra) with one layer per model.
              # You can visualize it directly:
              # plot(binary_maps)
              
  
  #//////////////  Project ensemble models (from single projections)   ///////////////////# 
  
  #EFpres <- BIOMOD_EnsembleForecasting( #bm.em = myBiomodEM, 
                                        # bm.proj = preds_pres, 
    #proj.name = 'CurrentEM',
    #models.chosen = 'all',
    #metric.binary = 'all',
    #metric.filter = 'all') #project ensemble to the present
  #EFpres_preds <- get_predictions(EFpres) - #don't need to do this step, as done below in output 
  #EFpres@proj.out@val
  
  #//////////////   Project ensemble models (building single projections)///////////////////# 
  EFpres <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                               proj.name = 'CurrentEM',
                                               new.env = as.data.frame(Clim),
                                               models.chosen = 'all',
                                               metric.binary = 'TSS',
                                               metric.filter = 'TSS')
  
  EFpres_preds <- get_predictions(EFpres) #don't need to do this step, as done below in output 
  EFpres@proj.out@val
  #unique(EFpres@proj.out@val$filtered.by)
  
  # Output the suitability values 
  filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_forecast_out/EFpres_out.csv"
  write_csv(EFpres@proj.out@val, file=filename, append = sp!=sp_names[1])
  
  binary_maps_cur <- get_predictions(EFpres, as.data.frame = FALSE, binary = TRUE)
  filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Binary_data/binary_data_current_ensemble.csv"
  write_csv(binary_maps_cur, file=filename, append = sp!=sp_names[1])
  #df <- as.data.frame(EFpres@proj.out@val)
  #unique(df$filtered.by)
 
  #table(df$filtered.by)
 

  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
  #  #  #  #  #  #  #  #  #  #  #  #  #  Projecting onto future climate (Single and ensemble )  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
  
  print("made it to future climate")

  scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
  # save them to individual CSVS, using the correct variables you chose. 
  for(v in scenarios) {
      
      names(Clim)
      Clim_fut=Clim # Using the Clim dataframe but changing to clim 
      
      #file_path <-read.csv("/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/ssp245_2041_2070.csv", row.names =1)
      
        file_path <-paste0("/Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/",v,"_2041_2070.csv") 
        Variables <- read.csv(file_path, row.names = 1)
        
        Variables <-  filter(Variables, row.names(Variables) %in% row.names(Clim))
      
   #row.names(ssp126_2041_2070) <- ssp126_2041_2070$ID1
   #unique(ssp126_2041_2070$ID1)
   #filename <- "//Users/catherinemercer/Downloads/ClimateAP/ClimateAP_v310/10_02_26/ssp126_2041_2070.csv"
   #write_csv(ssp126_2041_2070, file=filename, row.nam)
  
  
  # ssp245_2041_2070 <- GCMs %>% filter(GCM=="8GCMs_ensemble_ssp245_2041-2070.gcm")%>% 
  #  select(-GCM, -ID2)
  names(Clim_fut)
      Clim_fut$DD5 <- Variables$DD5
      Clim_fut$TD <- Variables$TD
      Clim_fut$MAP <- Variables$MAP
      Clim_fut$CMD <- Variables$CMD
 
      
      
      #//////////////  Single Models ///////////////////# 
      
      preds_fut <- BIOMOD_Projection(   bm.mod = models,
        new.env = Clim_fut,
       proj.name = 'CurrentRangeSize',
       selected.models = 'all',
        binary.meth = 'TSS'
      )
      
    
      #//////////////  Ensemble  ///////////////////# 

      EFfut <- BIOMOD_EnsembleForecasting( bm.em = myBiomodEM, 
                                           bm.proj = preds_fut,
                                           binary.meth ='TSS') #project ensemble to the future
            
      #EFpres_preds <- get_predictions(EFpres)
      fut_projection_df <- get_predictions(EFfut,  as.data.frame = TRUE)
      fut_projection_df$Species  <- sp
      fut_projection_df$Scenario <- v

    
      eval_scores <- get_evaluations(models)
      
      
      binary_maps_fut <- get_predictions(EFfut, as.data.frame = FALSE, binary = TRUE)
      binary_maps_fut$scenario <-v
      filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Binary_data/binary_data_future_ensemble.csv"
      write_csv(binary_maps_fut, file=filename, append = sp!=sp_names[1])

      # EFfut@models.out
      filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_forecast_out/fut_projection_df_2040_71.csv"
      write_csv(fut_projection_df, file=filename, append = sp!=sp_names[1])
      unique(fut_projection_df$Scenario)

      
      print("got to the binary maps and range change")
      binary_maps_fut
      binary_maps_cur
      
      threshold <- 500  # your TSS cutoff
      
      binary_maps_cur$current_bin <- ifelse(binary_maps_cur$pred>= threshold, 1, 0)
      binary_maps_fut$future_bin <- ifelse(binary_maps_fut$pred>= threshold, 1, 0)
      
      df_change <- NULL
      df_change$change[ binary_maps_cur$current_bin==1 & binary_maps_fut$future_bin==1] <- "stable"
      df_change$change[ binary_maps_cur$current_bin==1 & binary_maps_fut$future_bin==0] <- "loss"
      df_change$change[ binary_maps_cur$current_bin==0 & binary_maps_fut$future_bin==1] <- "gain"
      
      
     
    change <- as.data.frame(prop.table(table(df_change$change)))
    change$sp <- sp
    change$scenario <- NULL
    change$scenario <- v
    change <- change[, c("sp", "scenario", "Var1", "Freq")]
    
     table(df_change$change)
     change$current_range <- sum(binary_maps_cur$current_bin)
     change$future_range  <- sum(binary_maps_fut$future_bin)
     
     change$percent_change <- ((change$future_range - change$current_range ) / change$current_range ) * 100
    
     
    
    filename <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Binary_data/RangeChange.csv"
     write_csv(change, file=filename, append = sp!=sp_names[1])
     
      gc()
      print("End of loop reached for species....")
      print(sp)
      print("glm.fit: fitted probabilities numerically 0 or 1 occurre warning - likely small sample size and due to overfitting with 10 predictors, hopefully ensemble will remove overrconfident glm values - e.g. prediction of 1")
      
    
    }
  }







