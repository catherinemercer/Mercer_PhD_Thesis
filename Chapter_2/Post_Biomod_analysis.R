rm(list=ls()) # clear workspace

setwd("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Models")


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
library(tidyr)
library(tidyterra)
library(Hmisc)
library(lsr)
library(tidyverse)
  
  
  # https://ourcodingclub.github.io/tutorials/dataviz-beautification/
  theme_niwot <- function(){
    theme_bw() +
      theme(#text = element_text(family = "Arial"),
            axis.text = element_text(size = 16), 
            axis.title = element_text(size = 18),
            axis.line.x = element_line(color="black"), 
            axis.line.y = element_line(color="black"),
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),                                        panel.background = element_rect(colour='black'),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),  
            plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
            plot.title = element_text(size = 18, vjust = 1, hjust = 0),
            legend.text = element_text(size = 12),          
            legend.title = element_blank(),                              
            legend.position = c(0.95, 0.15), 
            legend.key = element_blank(),
            legend.background = element_rect(color = "black", 
                                             fill = "transparent", 
                                             size = 2, linetype = "blank"))
  }
  

#🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨

joined_demography_climate <- read_csv("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Post_SDM_Analysis/joined_demography_climate.csv")

names(joined_demography_climate)

#🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨 Accuracy measures 

names(joined_demography_climate)
boxplot_statistics<- joined_demography_climate %>% pivot_longer(KAPPA:TSS, names_to = "Accuracy_Metric", values_to = "Accuracy_value" )

boxplot_statistics_2 <- boxplot_statistics %>% pivot_longer(zones:Elev_class, names_to = "Class_type", values_to = "ClassValue" )  %>% 
  filter(!(ClassValue=="NA"))

boxplot_statistics_2$ClassValue <- recode(
  boxplot_statistics_2$ClassValue,
  "zonal" = "Zonal",
   "azonal"= "Azonal",
  "sub_montane" ="Sub-Montane",
  "montane" = "Montane",
   "high_montane"="Upper-Montane",
)

boxplot_statistics_2$ClassValue <- factor(boxplot_statistics_2$ClassValue, levels = c( "Sub-Montane", "Montane", "Upper-Montane", "Zonal", "Azonal",))

boxplot_statistics_2$Class_type <- recode(
  boxplot_statistics_2$Class_type,
  "Elev_class" = "Elevation Classes",
  "zones" = "Zonal Classes",
  
)

boxplot_statistics_2$ClassValue 

boxplot_statistics_2$Accuracy_Metric <- factor(boxplot_statistics_2$Accuracy_Metric, levels = c("TSS","ROC","KAPPA"))
  
bp <- ggplot(
  boxplot_statistics_2,
  aes(x = ClassValue, y = Accuracy_value, fill = ClassValue)
) +
  geom_boxplot(colour = "black") +
  facet_grid(
    ~Accuracy_Metric,
    scales = "free"
  ) +
  scale_fill_manual(
    values = c(
      "Sub-Montane"   = "#238b45",
      "Montane"       = "#74c476",
      "Upper-Montane" = "#c7e9c0",
      "Zonal"         = "#2171b5",
      "Azonal"        = "#a6bddb"
    )
  ) +
  labs(
    x = "Class",
    y = "Accuracy Score"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(
    legend.position = "none"
  ) + 
  theme(text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1))

bp


Summary statistics

#🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨 Variable importance


# Change in variable importance by elevation class and zonation + climate variable

# Think you just need to do a left join of the dataframe

VariableImportance <- read_csv("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Ensemble_Bmd_model_VarImport/EnsembleVariable_importance.csv")
print(VariableImportance)
unique(VariableImportance$algo)
nrow(VariableImportance)# 130880


VariableImportance <- VariableImportance  %>% group_by(sp,expl.var)  %>% 
  summarise(mean_var_imp=mean(var.imp))

VariableImportance$expl.var <- as.factor(VariableImportance$expl.var)
levels(VariableImportance$expl.var )



VariableImportance$expl.var <- factor(
  VariableImportance$expl.var,
  levels = c(
    "DD5",
    "CMD",
    "TD",
    "MAP",
    "eastness",
    "northness",
    "MeanCurvature",
    "slope",
    "TopographicPositionIndex",
    "TWI"
  ),
  labels = c(
    "Growing Degree Days",
    "Climate Moisture Deficit",
    "Temperature Seasonality",
    "Mean Annual Precipitation",
    "Eastness",
    "Northness",
    "Mean Curvature",
    "Slope",
    "Topographic Position Index",
    "Topographic Water Index"
  )
)


# 🌎🌎🌎🌎🌎㆞ making the plots

VariableImportance_Joined<- joined_demography_climate   %>% 
  left_join(VariableImportance, by="sp")


VarImport_boxplot_statistics_2 <- VariableImportance_Joined %>% pivot_longer(zones:Elev_class, names_to = "Class_type", values_to = "ClassValue" )  %>% 
  filter(!(ClassValue=="NA"))

VarImport_boxplot_statistics_2$ClassValue <- recode(
  VarImport_boxplot_statistics_2$ClassValue,
  "zonal" = "Zonal",
  "azonal"= "Azonal",
  "sub_montane" ="Sub-Montane",
  "montane" = "Montane",
  "high_montane"="Upper-Montane",
)

VarImport_boxplot_statistics_2$ClassValue <- as.factor(VarImport_boxplot_statistics_2$ClassValue)
VarImport_boxplot_statistics_2$ClassValue <- factor(VarImport_boxplot_statistics_2$ClassValue, levels = c("Sub-Montane","Montane","Upper-Montane", "Zonal", "Azonal"))


table(VarImport_boxplot_statistics_2$Class)
#ALL class types
bp <- ggplot(
  VarImport_boxplot_statistics_2,
  aes(x = ClassValue, y = mean_var_imp, fill = ClassValue)
) +
  geom_boxplot(colour = "black") +
  facet_wrap(
    ~expl.var,
    ncol = 5,
    nrow = 2,
    scales = "free"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(
    values = c(
      "Sub-Montane"   = "#238b45",
      "Montane"       = "#74c476",
      "Upper-Montane" = "#c7e9c0",
      "Zonal"         = "#2171b5",
      "Azonal"        = "#a6bddb"
    )
  ) +
  labs(
    x = "Class",
    y = "Variable Importance"
  )  +
  theme_bw()+
  theme(text = element_text(size = 18),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))


#🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨 1981-2010 vs adult

joined_demography_climate <- read_csv("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/1 - Project Folder Template/3-ExperimentOne/4-Outputs/SDM_Outputs/110225/Post_SDM_Analysis/joined_demography_climate.csv")

names(joined_demography_climate)


  
#🟧🟨 LOWER
hist(joined_demography_climate$lower_elev_adult)
hist(joined_demography_climate$l_current)

hist(joined_demography_climate$lower_elev_adult,
     col = rgb(1, 0, 0, 0.5),
     xlim = range(c(joined_demography_climate$lower_elev_adult,
                    joined_demography_climate$l_current)),
     main = "",
     xlab = "Lower elevation (m)")

hist(joined_demography_climate$l_current,
     col = rgb(0, 0, 1, 0.5),
     add = TRUE)

# RUN FROM HERE TO PATCHWORK PLOTTING
# LOWER
set.seed(123)
mod <- lm(lower_elev_adult ~ l_current, data = joined_demography_climate) 
lower_current_mod <- mod
summary(lower_current_mod) #🟩🟩🟩🟩
smaller <- joined_demography_climate$l_current < 
  joined_demography_climate$lower_elev_adult

mod_df <- data.frame(
  lower_elev_adult = joined_demography_climate$lower_elev_adult,
  l_current = joined_demography_climate$l_current,
  smaller = smaller
)

Lower_adult_current_plot <- ggplot(
  mod_df, 
  aes(x = l_current, y = lower_elev_adult)
) +
  geom_point(aes(colour = smaller), alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "black",
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#54278f", `TRUE` = "#54278f"),
    guide = "none"
  ) +
  labs(
    x = "SDM predicted elevation limit (m)",
    y = "Elevation limit of adult trees (m)"
  ) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 4000)) +
  coord_equal() +
  xlim(0,3000)+
  ylim(0,3000)+
  theme_classic() +
  theme(text = element_text(size = 20))



# MEDIAN
set.seed(123)
mod <- lm(med_elev_adult ~ m_current, data = joined_demography_climate) 
median_current_mod <- mod
summary(median_current_mod)#🟩🟩🟩🟩
smaller <- joined_demography_climate$m_current < 
  joined_demography_climate$med_elev_adult

mod_df <- data.frame(
  med_elev_adult = joined_demography_climate$med_elev_adult,
  m_current = joined_demography_climate$m_current,
  smaller = smaller
)

median_adult_current_plot <- ggplot(
  mod_df, 
  aes(x = m_current, y = med_elev_adult)
) +
  geom_point(aes(colour = smaller), alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "black",
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#807dba", `TRUE` = "#807dba"),
    guide = "none"
  ) +
  labs(
    x = "SDM predicted elevation limit (m)",
    y = ""
  ) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 4000)) +
  coord_equal() +
  xlim(0,3000)+
  ylim(0,3000)+
  theme_classic() +
  theme(text = element_text(size = 20))

# UPPER
set.seed(123)
mod <- lm(upper_elev_adult ~ u_current, data = joined_demography_climate) 
upper_current_mod <- mod
summary(upper_current_mod) #🟩🟩🟩🟩

smaller <- joined_demography_climate$u_current < 
  joined_demography_climate$upper_elev_adult

mod_df <- data.frame(
  upper_elev_adult = joined_demography_climate$upper_elev_adult,
  u_current = joined_demography_climate$u_current,
  smaller = smaller
)

upper_adult_current_plot <- ggplot(
  mod_df, 
  aes(x = u_current, y = upper_elev_adult)
) +
  geom_point(aes(colour = smaller), alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "black",
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#bcbddc", `TRUE` = "#bcbddc"),
    guide = "none"
  ) +
  labs(
    x = "SDM predicted elevation limit (m)",
    y = ""
  ) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 4000)) +
  coord_equal() +
  theme_classic() +
  xlim(0,3000)+
  ylim(0,3000)+
  theme(text = element_text(size = 20))

library(patchwork)
Lower_adult_current_plot + median_adult_current_plot + upper_adult_current_plot

combined_plot <- (
  Lower_adult_current_plot +
    median_adult_current_plot +
    upper_adult_current_plot 
) +
  plot_layout(ncol = 3)

combined_plot




#🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨 2041-2070 SDM predictions vs Juveniles 
raw_data <- joined_demography_climate %>% 
  select(sp, zones:upper_elev_juvenile) %>% 
  select(-m_current, -u_current, -l_current, -med_elev_adult,-lower_elev_adult, -upper_elev_adult) %>% 
  rename(m_126 = "m_f204170_126",
         m_245= "m_f204170_245",
         m_370= "m_f204170_370", 
         m_585= "m_f204170_585",
         u_126 ="u_f204170_126", 
         u_245 ="u_f204170_245",
         u_370 = "u_f204170_370", 
         u_585 = "u_f204170_585", 
         l_126 = "l_f204170_126", 
         l_245 = "l_f204170_245", 
         l_370 = "l_f204170_370", 
         l_585 = "l_f204170_585",
         l_juv = "lower_elev_juvenile",
         m_juv = "med_elev_juvenile", 
         u_juv = "upper_elev_juvenile"
  )


#🟧🟨 


LOWER_raw_data <- raw_data %>%  
  select(sp, l_126, l_245, l_370, l_585, l_juv) %>%  
  pivot_longer(
    cols = c(l_126, l_245, l_370, l_585),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      l_126 = "126",
      l_245 = "245",
      l_370 = "370",
      l_585 = "585"
    ),
    smaller = l_juv < Value
  )
#select(sp,l_126,l_245, l_370, l_585, l_juv) %>%  
# pivot_longer(l_126:l_585, names_to = "Scenario", values_to = "Value") 

LOWER_raw_data[LOWER_raw_data=="l_126"]<-"126"
LOWER_raw_data[LOWER_raw_data=="l_245"]<-"245"
LOWER_raw_data[LOWER_raw_data=="l_370"]<-"370"
LOWER_raw_data[LOWER_raw_data=="l_585"]<-"585"

names(LOWER_raw_data)

set.seed(123)
lower_models <- LOWER_raw_data %>%
  group_by(Scenario) %>%
  group_modify(~ {
    mod <- lm(l_juv ~ Value, data = .x)
    print(summary(mod))
    s <- summary(mod)
    
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()

lower_models

# LOWER
p_l <- ggplot(LOWER_raw_data, aes(x = Value, y = l_juv)) +
  geom_point(aes(colour = smaller), alpha = 0.7, size = 2) +
  
  geom_smooth(
    aes(
      group = Scenario,
      colour = Scenario
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.15
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = c(
      "126" = "#fecc5c",
      "245" = "#fd8d3c",
      "370" = "#f03b20",
      "585" = "#bd0026",
      "TRUE" = "#54278f",
      "FALSE" = "#54278f"
      
      #fecc5c
      #fd8d3c
      #f03b20
      #bd0026
    ),
    guide = "none"
  ) +
  
  labs(
    x = "SDM Predicted Elevation Limit (m)",
    y = "Elevation limit of juvenile trees (m)"
  ) +
  
coord_cartesian(xlim = c(0, 4000), ylim = c(0, 4000)) +
  coord_equal() +
  theme_classic() +
  xlim(0,3000)+
  ylim(0,3000)+
  theme(text = element_text(size = 20))


#🟧🟨 Median 


MEDIAN_raw_data <- raw_data %>%  
  select(sp, m_126, m_245, m_370, m_585, m_juv) %>%  
  pivot_longer(
    cols = c(m_126, m_245, m_370, m_585),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      m_126 = "126",
      m_245 = "245",
      m_370 = "370",
      m_585 = "585"
    ),
    smaller = m_juv < Value
  )

set.seed(123)
median_models <- MEDIAN_raw_data %>%
  group_by(Scenario) %>%
  group_modify(~ {
    mod <- lm(m_juv ~ Value, data = .x)
    s <- summary(mod)
    print(summary(mod))
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()
median_models

p_m <- ggplot(MEDIAN_raw_data, aes(x = Value, y = m_juv)) +
  geom_point(aes(colour = smaller), alpha = 0.7, size = 2) +
  
  geom_smooth(
    aes(
      group = Scenario,
      colour = Scenario
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.15
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = c(
      "126" = "#fecc5c",
      "245" = "#fd8d3c",
      "370" = "#f03b20",
      "585" = "#bd0026",
      "TRUE" = "#807dba",
      "FALSE" = "#807dba"
      
      #fecc5c
      #fd8d3c
      #f03b20
      #bd0026
    ),
    guide = "none"
  ) +
  
  labs(
    x = "SDM Predicted Elevation Limit (m)",
    y = ""
  ) +
  coord_cartesian(xlim = c(0, 4000), ylim = c(0, 4000)) +
  coord_equal() +
  theme_classic() +
  xlim(0,3000)+
  ylim(0,3000)+
  theme(text = element_text(size = 20))


#🟧🟨 Upper

UPPER_raw_data <- raw_data %>%  
  select(sp, u_126, u_245, u_370, u_585, u_juv) %>%  
  pivot_longer(
    cols = c(u_126, u_245, u_370, u_585),
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      u_126 = "126",
      u_245 = "245",
      u_370 = "370",
      u_585 = "585"
    ),
    smaller = u_juv < Value
  )

set.seed(123)
upper_models <- UPPER_raw_data %>%
  group_by(Scenario) %>%
  group_modify(~ {
    mod <- lm(u_juv ~ Value, data = .x)
    s <- summary(mod)
    print(summary(mod))
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()

upper_models
# UPPER

p_u <- ggplot(UPPER_raw_data, aes(x = Value, y = u_juv)) +
  geom_point(aes(colour = smaller), alpha = 0.7, size = 2) +
  
  geom_smooth(
    aes(
      group = Scenario,
      colour = Scenario
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.15
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = c(
      "126" = "#fecc5c",
      "245" = "#fd8d3c",
      "370" = "#f03b20",
      "585" = "#bd0026",
      "TRUE" = "#bcbddc",
      "FALSE" = "#bcbddc"
      
      #fecc5c
      #fd8d3c
      #f03b20
      #bd0026
    ),
    guide = "none"
  ) +
  labs(
    x = "SDM Predicted Elevation Limit (m)",
    y = ""
  ) +
coord_cartesian(xlim = c(0, 4000), ylim = c(0, 4000)) +
  coord_equal() +
  theme_classic() +
  xlim(0,3000)+
  ylim(0,3000)+
  theme(text = element_text(size = 20))

# All three plots
library(patchwork)
p_l + p_m + p_u



#Elevation Difference 🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨


joined_demography_climate_summary <- joined_demography_climate %>% 
  mutate(
    obs_change_median = med_elev_juvenile - med_elev_adult,
    
    pred_change_median_126 = m_f204170_126 - m_current, 
    pred_change_median_245 = m_f204170_245 - m_current,
    pred_change_median_370 = m_f204170_370 - m_current, 
    pred_change_median_585 = m_f204170_585 - m_current, 
    
    obs_change_upper = upper_elev_juvenile - upper_elev_adult, 
    
    pred_change_upper_126 = u_f204170_126 - u_current, 
    pred_change_upper_245 = u_f204170_245 - u_current, 
    pred_change_upper_370 = u_f204170_370 - u_current, 
    pred_change_upper_585 = u_f204170_585 - u_current, 
    
    obs_change_lower = lower_elev_juvenile - lower_elev_adult,
    
    pred_change_lower_126 = l_f204170_126 - l_current, 
    pred_change_lower_245 = l_f204170_245 - l_current,
    pred_change_lower_370 = l_f204170_370 - l_current, 
    pred_change_lower_585 = l_f204170_585 - l_current
  ) %>% 
  rename(
    m_c_126 = pred_change_median_126,
    m_c_245 = pred_change_median_245,
    m_c_370 = pred_change_median_370, 
    m_c_585 = pred_change_median_585,
    
    u_c_126 = pred_change_upper_126, 
    u_c_245 = pred_change_upper_245, 
    u_c_370 = pred_change_upper_370, 
    u_c_585 = pred_change_upper_585, 
    
    l_c_126 = pred_change_lower_126, 
    l_c_245 = pred_change_lower_245, 
    l_c_370 = pred_change_lower_370, 
    l_c_585 = pred_change_lower_585,
    
    l_c_juv = obs_change_lower,
    m_c_juv = obs_change_median, 
    u_c_juv = obs_change_upper
  )




LOWER_change <- joined_demography_climate_summary %>%  
  select(sp, l_c_126, l_c_245, l_c_370, l_c_585, l_c_juv) %>%  
  pivot_longer(
    l_c_126:l_c_585,
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      l_c_126 = "126",
      l_c_245 = "245",
      l_c_370 = "370",
      l_c_585 = "585"
    )
  ) %>%
  filter(sp != "Rhododendron.leptosanthum")


lower_change_models <- LOWER_change %>%
  group_by(Scenario) %>%
  group_modify(~ {
    set.seed(123)
    mod <- lm(l_c_juv ~ Value, data = .x)
    s <- summary(mod)
    print(summary(mod))
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()

lower_change_models




LOWER_change <- LOWER_change %>%
  mutate(
    above = l_c_juv > Value
  )

p_l_change <- ggplot(
  LOWER_change,
  aes(x = Value, y = l_c_juv)
) +
  
  # Points: above/below 1:1
  geom_point(
    aes(colour = above),
    alpha = 0.7,
    size = 2
  ) +
  
  # Regression line for each scenario
  geom_smooth(
    aes(
      colour = Scenario,
      group = Scenario
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.15
  ) +
  
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = c(
      "126" = "#fecc5c",
      "245" = "#fd8d3c",
      "370" = "#f03b20",
      "585" = "#bd0026",
      "FALSE" = "#54278f",
      "TRUE" = "#54278f"
    ),
    guide = "none"
  ) +
  
  labs(
    x = "SDM temporal elevation shift (m)",
    y = "Demographic elevation shift (m)"
  ) +
  ylim(-500, 500) +
  xlim(-500, 500) +
  coord_equal() +
  theme_classic() +
  theme(text = element_text(size = 20))




#🟧🟨 Median
MEDIAN_change <- joined_demography_climate_summary %>%  
  select(sp, m_c_126, m_c_245, m_c_370, m_c_585, m_c_juv) %>%  
  pivot_longer(
    m_c_126:m_c_585,
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      m_c_126 = "126",
      m_c_245 = "245",
      m_c_370 = "370",
      m_c_585 = "585"
    )
  ) %>%
  filter(sp != "Rhododendron.leptosanthum")


median_change_models <- MEDIAN_change %>%
  group_by(Scenario) %>%
  group_modify(~ {
    
    mod <- lm(m_c_juv ~ Value, data = .x)
    s <- summary(mod)
    print(s)
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()

median_change_models

MEDIAN_change <- joined_demography_climate_summary %>%  
  select(sp, m_c_126, m_c_245, m_c_370, m_c_585, m_c_juv) %>%  
  pivot_longer(
    m_c_126:m_c_585,
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      m_c_126 = "126",
      m_c_245 = "245",
      m_c_370 = "370",
      m_c_585 = "585"
    )
  ) %>%
  filter(sp != "Rhododendron.leptosanthum")


median_change_models <- MEDIAN_change %>%
  group_by(Scenario) %>%
  group_modify(~ {
    
    mod <- lm(m_c_juv ~ Value, data = .x)
    s <- summary(mod)
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()

median_change_models
MEDIAN_change <- MEDIAN_change %>%
  mutate(
    above = m_c_juv > Value
  )

p_m_change <- ggplot(
  MEDIAN_change,
  aes(x = Value, y = m_c_juv)
) +
  
  # Points: above/below 1:1
  geom_point(
    aes(colour = above),
    alpha = 0.7,
    size = 2
  ) +
  
  # Regression line for each scenario
  geom_smooth(
    aes(
      colour = Scenario,
      group = Scenario
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.15
  ) +
  
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = c(
      "126" = "#fecc5c",
      "245" = "#fd8d3c",
      "370" = "#f03b20",
      "585" = "#bd0026",
      "FALSE" = "#807dba",
      "TRUE" = "#807dba"
    ),
    guide = "none"
  ) +

  labs(
    x = "SDM temporal elevation shift (m)",
    y = ""
  ) +
  coord_equal() +
  ylim(-500, 500) +
  xlim(-500, 500) +
  coord_equal() +
  theme_classic() +
  theme(text = element_text(size = 20))


#🟧🟨 Upper
UPPER_change <- joined_demography_climate_summary %>%  
  select(sp, u_c_126, u_c_245, u_c_370, u_c_585, u_c_juv) %>%  
  pivot_longer(
    u_c_126:u_c_585,
    names_to = "Scenario",
    values_to = "Value"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      u_c_126 = "126",
      u_c_245 = "245",
      u_c_370 = "370",
      u_c_585 = "585"
    )
  ) %>%
  filter(sp != "Rhododendron.leptosanthum")


upper_change_models <- UPPER_change %>%
  group_by(Scenario) %>%
  group_modify(~ {
    
    mod <- lm(u_c_juv ~ Value, data = .x)
    s <- summary(mod)
    print(s)
    p <- coef(s)[2, 4]
    
    data.frame(
      intercept = coef(mod)[1],
      slope = coef(mod)[2],
      R2 = s$r.squared,
      Adj_R2 = s$adj.r.squared,
      p_value = p,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE      ~ "ns"
      )
    )
  }) %>%
  ungroup()

upper_change_models

UPPER_change <- UPPER_change %>%
  mutate(
    above = u_c_juv > Value
  )

p_u_change <- ggplot(
  UPPER_change,
  aes(x = Value, y = u_c_juv)
) +
  
  # Points: above/below 1:1
  geom_point(
    aes(colour = above),
    alpha = 0.7,
    size = 2
  ) +
  
  # Regression line for each scenario
  geom_smooth(
    aes(
      colour = Scenario,
      group = Scenario
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.15
  ) +
  
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = c(
      "126" = "#fecc5c",
      "245" = "#fd8d3c",
      "370" = "#f03b20",
      "585" = "#bd0026",
      "FALSE" = "#bcbddc",
      "TRUE" = "#bcbddc"
    ),
    guide = "none"
  ) +
  
  labs(
    x = "SDM temporal elevation shift (m)",
    y = ""
  ) +
  coord_equal() +
  ylim(-500, 500) +
  xlim(-500, 500) +
  coord_equal() +
  theme_classic() +
  theme(text = element_text(size = 20))



p_l_change + p_m_change +p_u_change



#🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧🟨  Barplots Change over time 

# I think i used the same code for zonal/elevation groupings and just amended the code

joined_demography_climate_summary <- joined_demography_climate %>% 
  mutate(
    obs_change_median = med_elev_juvenile - med_elev_adult,
    
    pred_change_median_126 = m_f204170_126 - m_current, 
    pred_change_median_245 = m_f204170_245 - m_current,
    pred_change_median_370 = m_f204170_370 - m_current, 
    pred_change_median_585 = m_f204170_585 - m_current, 
    
    obs_change_upper = upper_elev_juvenile - upper_elev_adult, 
    
    pred_change_upper_126 = u_f204170_126 - u_current, 
    pred_change_upper_245 = u_f204170_245 - u_current, 
    pred_change_upper_370 = u_f204170_370 - u_current, 
    pred_change_upper_585 = u_f204170_585 - u_current, 
    
    obs_change_lower = lower_elev_juvenile - lower_elev_adult,
    
    pred_change_lower_126 = l_f204170_126 - l_current, 
    pred_change_lower_245 = l_f204170_245 - l_current,
    pred_change_lower_370 = l_f204170_370 - l_current, 
    pred_change_lower_585 = l_f204170_585 - l_current
  ) %>% 
  rename(
    m_c_126 = pred_change_median_126,
    m_c_245 = pred_change_median_245,
    m_c_370 = pred_change_median_370, 
    m_c_585 = pred_change_median_585,
    
    u_c_126 = pred_change_upper_126, 
    u_c_245 = pred_change_upper_245, 
    u_c_370 = pred_change_upper_370, 
    u_c_585 = pred_change_upper_585, 
    
    l_c_126 = pred_change_lower_126, 
    l_c_245 = pred_change_lower_245, 
    l_c_370 = pred_change_lower_370, 
    l_c_585 = pred_change_lower_585,
    
    l_c_juv = obs_change_lower,
    m_c_juv = obs_change_median, 
    u_c_juv = obs_change_upper
  )

library(tidyverse)
joined_demography_climate$zones
 


library(tidyverse)
library(tidyverse)

# 1. Clean, filter, and pivot the data
raw_data_adapted <- joined_demography_climate %>%
  filter(sp != "Rhododendron.leptosanthum") %>%
  filter(Elev_class %in% c("sub_montane", "montane", "high_montane")) %>%
  
  mutate(
    l_demographic = lower_elev_juvenile - lower_elev_adult,
    m_demographic = med_elev_juvenile - med_elev_adult,
    u_demographic = upper_elev_juvenile - upper_elev_adult,
    
    l_126 = l_f204170_126 - l_current, l_245 = l_f204170_245 - l_current,
    l_370 = l_f204170_370 - l_current, l_585 = l_f204170_585 - l_current,
    
    m_126 = m_f204170_126 - m_current, m_245 = m_f204170_245 - m_current,
    m_370 = m_f204170_370 - m_current, m_585 = m_f204170_585 - m_current,
    
    u_126 = u_f204170_126 - u_current, u_245 = u_f204170_245 - u_current,
    u_370 = u_f204170_370 - u_current, u_585 = u_f204170_585 - u_current
  ) %>%
  select(sp, Elev_class, matches("^([lmu])_(demographic|126|245|370|585)$")) %>%
  
  pivot_longer(
    cols = -c(sp, Elev_class),
    names_to = c("Limit", "Group"),
    names_pattern = "([mul])_(.*)",
    values_to = "Change"
  ) %>%
  mutate(
    Limit = recode(Limit, l = "Lower", m = "Median", u = "Upper"),
    Limit = factor(Limit, levels = c("Lower", "Median", "Upper")),
    
    Group = factor(
      Group,
      levels = c("demographic", "126", "245", "370", "585"),
      # CHANGED: "Observed Demography" label shortened to "Demography"
      labels = c("Demography", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    ),
    
    # Matching the updated Demography factor label
    Source = if_else(Group == "Demography", "Demographic", "SDM Prediction"),
    
    Elev_class = recode(
      Elev_class,
      sub_montane  = "Sub-Montane",
      montane      = "Montane",
      high_montane = "Upper-Montane"
    ),
    Elev_class = factor(Elev_class, levels = c("Sub-Montane", "Montane", "Upper-Montane"))
  )


# 2. Compute summaries (Mean and Standard Error) 🟧🟧🟨🟨

summary_data <- raw_data_adapted %>%
  group_by(Elev_class, Limit, Group, Source) %>%
  summarise(
    Mean_Change = mean(Change, na.rm = TRUE),
    SE = sd(Change, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
citation("mgcv")



# 3. Create the side-by-side bar chart
p_limit_change_adapted <- ggplot(
  summary_data,
  aes(x = Group, y = Mean_Change, fill = Source)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(ymin = Mean_Change - SE, ymax = Mean_Change + SE),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "gray40"
  ) +
  
  facet_grid(
    Limit ~ Elev_class
  ) +
  
  coord_cartesian(ylim = c(-150, 350)) +
  
  scale_fill_manual(
    values = c(
      "Demographic"    = "#d95f02",
      "SDM Prediction" = "#41ab5d"
    )
  ) +
  labs(
    x = "Demographic and SDM predicted change",
    y = "Elevation Limit Change (m)",
    fill = "Prediction Source"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_limit_change_adapted



# Tukey test 🟧🟧🟨🟨

library(emmeans)
library(multcomp)
library(multcompView)

# 1. Fit the linear model using your RAW dataset (not the summary data frame)
set.seed(123)
model <- lm(Change ~ Elev_class * Limit * Group, data = raw_data_adapted)


summary(model)
# 2. Extract Estimated Marginal Means grouped by your facets
emm_results <- emmeans(model, specs = ~ Group | Elev_class * Limit)

# 3. Generate the Compact Letter Display (Tukey-adjusted)
cld_letters <- cld(emm_results, adjust = "Tukey", Letters = letters, alpha = 0.05)

# 4. Clean up the output to merge with your summary data
# Note: cld() outputs a column named ".group" containing the letters
letters_df <- as.data.frame(cld_letters) %>%
  mutate(.group = trimws(.group)) %>% # Strip trailing/leading spaces
  dplyr::select(Elev_class, Limit, Group, letter = .group)

summary_data <- summary_data %>%
  left_join(letters_df, by = c("Elev_class", "Limit", "Group"))


# Place this at the end of your ggplot chain, right before theme_bw()
  
  
  
  # 3. Create the side-by-side bar chart
  p_limit_change_adapted <- ggplot(
    summary_data,
    aes(x = Group, y = Mean_Change, fill = Source)
  ) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(ymin = Mean_Change - SE, ymax = Mean_Change + SE),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "gray40"
  ) +
  
  facet_grid(
    Limit ~ Elev_class
  ) +
  
  coord_cartesian(ylim = c(-150, 350)) +
  
  scale_fill_manual(
    values = c(
      "Demographic"    = "#d95f02",
      "SDM Prediction" = "#41ab5d"
    )
  ) +
  labs(
    x = "Data Source",
    y = "Elevation Limit Change (m)",
    fill = ""
  ) +
  theme_bw() +
    geom_text(
      data = summary_data,
      aes(
        x = Group,
        # Places letters slightly above the positive SE bar or below negative SE bar
        y = if_else(Mean_Change >= 0, Mean_Change + SE + 15, Mean_Change - SE - 25), 
        label = letter
      ),
      position = position_dodge(width = 0.8),
      vjust = if_else(summary_data$Mean_Change >= 0, 0, 1),
      size = 5,
      fontface = "bold"
    ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(), 
    legend.text=element_text(size=14)
  )

p_limit_change_adapted


#####🟧🟧🟨🟨🟧🟧🟨🟨🟧🟧🟨🟨🟧🟧🟨🟨🟧🟧🟨🟨🟧🟧🟨🟨
# 1. Clean, filter, and pivot the data
names(raw_data_adapted)
joined_demography_climate$zones
names(joined_demography_climate)
raw_data_adapted <- joined_demography_climate %>%
  filter(sp != "Rhododendron.leptosanthum") %>%
  filter(zones %in% c("zonal", "azonal")) %>%
  
  mutate(
    l_demographic = lower_elev_juvenile - lower_elev_adult,
    m_demographic = med_elev_juvenile - med_elev_adult,
    u_demographic = upper_elev_juvenile - upper_elev_adult,
    
    l_126 = l_f204170_126 - l_current, l_245 = l_f204170_245 - l_current,
    l_370 = l_f204170_370 - l_current, l_585 = l_f204170_585 - l_current,
    
    m_126 = m_f204170_126 - m_current, m_245 = m_f204170_245 - m_current,
    m_370 = m_f204170_370 - m_current, m_585 = m_f204170_585 - m_current,
    
    u_126 = u_f204170_126 - u_current, u_245 = u_f204170_245 - u_current,
    u_370 = u_f204170_370 - u_current, u_585 = u_f204170_585 - u_current
  ) %>%
  select(sp, zones, matches("^([lmu])_(demographic|126|245|370|585)$")) %>%
  
  pivot_longer(
    cols = -c(sp, zones),
    names_to = c("Limit", "Group"),
    names_pattern = "([mul])_(.*)",
    values_to = "Change"
  ) %>%
  mutate(
    Limit = recode(Limit, l = "Lower", m = "Median", u = "Upper"),
    Limit = factor(Limit, levels = c("Lower", "Median", "Upper")),
    
    Group = factor(
      Group,
      levels = c("demographic", "126", "245", "370", "585"),
      # CHANGED: "Observed Demography" label shortened to "Demography"
      labels = c("Demography", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    ),
    
    # Matching the updated Demography factor label
    Source = if_else(Group == "Demography", "Demographic", "SDM Prediction"),
    
    zones = recode(
      zones,
      zonal  = "Zonal",
      azonal = "Azonal"
    ),
    zones = factor(zones, levels = c("Zonal", "Azonal"))
  )


# 2. Compute summaries (Mean and Standard Error) 🟧🟧🟨🟨

summary_data <- raw_data_adapted %>%
  group_by(zones, Limit, Group, Source) %>%
  summarise(
    Mean_Change = mean(Change, na.rm = TRUE),
    SE = sd(Change, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
 citation("mgcv")

print(summary_data)

# 3. Create the side-by-side bar chart
p_limit_change_adapted <- ggplot(
  summary_data,
  aes(x = Group, y = Mean_Change, fill = Source)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(ymin = Mean_Change - SE, ymax = Mean_Change + SE),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "gray40"
  ) +
  
  facet_grid(
    Limit ~ zones
  ) +
  
  coord_cartesian(ylim = c(-150, 350)) +
  
  scale_fill_manual(
    values = c(
      "Demographic"    = "#fa9fb5",
      "SDM Prediction" = "#1c9099"
    )
  ) +
  labs(
    x = "Demographic and SDM predicted change",
    y = "Elevation Limit Change (m)",
    fill = "Prediction Source"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_limit_change_adapted



# Tukey test 🟧🟧🟨🟨

library(emmeans)
library(multcomp)
library(multcompView)

# 1. Fit the linear model using your RAW dataset (not the summary data frame)
set.seed(123)
model <- lm(Change ~ zones * Limit * Group, data = raw_data_adapted)


summary(model)
# 2. Extract Estimated Marginal Means grouped by your facets
emm_results <- emmeans(model, specs = ~ Group | zones * Limit)

# 3. Generate the Compact Letter Display (Tukey-adjusted)
cld_letters <- cld(emm_results, adjust = "Tukey", Letters = letters, alpha = 0.05)

# 4. Clean up the output to merge with your summary data
# Note: cld() outputs a column named ".group" containing the letters
letters_df <- as.data.frame(cld_letters) %>%
  mutate(.group = trimws(.group)) %>% # Strip trailing/leading spaces
  dplyr::select(zones, Limit, Group, letter = .group)

summary_data <- summary_data %>%
  left_join(letters_df, by = c("zones", "Limit", "Group"))

# Place this at the end of your ggplot chain, right before theme_bw()

# 3. Create the side-by-side bar chart
p_limit_change_adapted <- ggplot(
  summary_data,
  aes(x = Group, y = Mean_Change, fill = Source)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(ymin = Mean_Change - SE, ymax = Mean_Change + SE),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "gray40"
  ) +
  
  facet_grid(
    Limit ~ zones
  ) +
  
  coord_cartesian(ylim = c(-150, 350)) +
  
  scale_fill_manual(
    values = c(
      "Demographic"    = "#fa9fb5",
      "SDM Prediction" = "#1c9099"
    )
  ) +
  labs(
    x = "Data Source",
    y = "Elevation Limit Change (m)",
    fill = ""
  ) +
  theme_bw() +
  geom_text(
    data = summary_data,
    aes(
      x = Group,
      # Places letters slightly above the positive SE bar or below negative SE bar
      y = if_else(Mean_Change >= 0, Mean_Change + SE + 15, Mean_Change - SE - 25), 
      label = letter
    ),
    position = position_dodge(width = 0.8),
    vjust = if_else(summary_data$Mean_Change >= 0, 0, 1),
    size = 5,
    fontface = "bold"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(), 
    legend.text=element_text(size=14)
  )

p_limit_change_adapted
