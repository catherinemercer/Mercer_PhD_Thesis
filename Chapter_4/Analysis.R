rm(list=ls())
library(terra)
library(ggplot2)
library(patchwork)
library(ggridges)
library(leaflet)
library(tidyr)
library(see)
library(dplyr)
library(mgcv)
library(gratia)
library(tidyverse)
print("You got this dude")

#  🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Forest mask
#LandcoverMap <- rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/January_25/12_01_25/2022_only/Landsat_2002_L5/Lansat_2002.tif")

# Forest_class <- LandcoverMap$Class==0
# Forest_class_TRUE <- subst(Forest_class, FALSE, NA)
#plet(Forest_class_TRUE, y=1, alpha=0.8, 
# tiles=c("Esri.WorldImagery"))


# Correct function for reading shapefiles (.shp)
Landcover_vect <- vect("/Users/catherinemercer/TaiwanForest_VegetationMap.shp")

# 2. Extract only the rows where the class_en column equals "Forest"
Forest_class_TRUE_vect <- Landcover_vect[Landcover_vect$class_en == "Forest", ]
names(Landcover_vect)
plot(Forest_class_TRUE_vect)


#🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Filter by year

LandtrendR <- rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/LandtrendR/lt-loss_map_060325.tif")
crs(LandtrendR)
hist(LandtrendR$mag) # skewed
hist(log1p(LandtrendR$mag)) # more normal

LandtrendR <- mask(LandtrendR, Forest_class_TRUE_vect)

plet(LandtrendR,tiles=c("Esri.WorldImagery"))

#🟧🟨🟧🟨🟧🟨🟧🟨 Define change and no change pixels 
# 🔴🔴🔴🔴 Change
set.seed(123)
change_20132017 <- LandtrendR$yod >= 2013 & LandtrendR$yod <= 2017
plet(change_20132017,tiles=c("Esri.WorldImagery"))
# 🟡🟡🟡🟡 No Change
nochange_20132017 <- ifel(change_20132017 == 0, 1, NA)
plet(nochange_20132017,tiles=c("Esri.WorldImagery"))

# 🔴🔴🔴🔴 Change
set.seed(123)
change_20182022 <- LandtrendR$yod >= 2018 & LandtrendR$yod <= 2022
plet(change_20182022,tiles=c("Esri.WorldImagery"))
# 🟡🟡🟡🟡 No Change
nochange_20182022 <- ifel(change_20182022 == 0, 1, NA)
plet(nochange_20182022,tiles=c("Esri.WorldImagery"))

  #Recalculate pixels 

# Count changed pixels in each period
change_1 <- global(change_20132017, "sum", na.rm = TRUE)[1,1]
change_2 <- global(change_20182022, "sum", na.rm = TRUE)[1,1]

# Count no-change pixels
nochange_1 <- global(nochange_20132017, "sum", na.rm = TRUE)[1,1]
nochange_2 <- global(nochange_20182022, "sum", na.rm = TRUE)[1,1]

# Contingency table
chi_table <- matrix(
  c(change_1, nochange_1,
    change_2, nochange_2),
  nrow = 2,
  byrow = TRUE
)


rownames(chi_table) <- c("2013-2017", "2018-2022")
colnames(chi_table) <- c("Change", "No Change")

chi_result <- chisq.test(chi_table)
chi_result

percentage_change <- prop.table(chi_table, margin = 1) * 100
percentage_change


library(terra)

# Calculate area of each pixel in hectares
pixel_area <- cellSize(LandtrendR, unit = "ha")
# Period 1: 2013–2017
change_area_1 <- global(
  pixel_area * change_20132017,
  "sum", na.rm = TRUE
)[1,1]

nochange_area_1 <- global(
  pixel_area * nochange_20132017,
  "sum", na.rm = TRUE
)[1,1]

# Period 2: 2018–2022
change_area_2 <- global(
  pixel_area * change_20182022,
  "sum", na.rm = TRUE
)[1,1]

nochange_area_2 <- global(
  pixel_area * nochange_20182022,
  "sum", na.rm = TRUE
)[1,1]


area_summary <- data.frame(
  Period = c("2013-2017", "2018-2022"),
  Change_ha = c(change_area_1, change_area_2),
  NoChange_ha = c(nochange_area_1, nochange_area_2)
)

area_summary$Total_ha <- 
  area_summary$Change_ha + area_summary$NoChange_ha

area_summary$Change_km2 <- area_summary$Change_ha / 100
area_summary$NoChange_km2 <- area_summary$NoChange_ha / 100
area_summary$Total_km2 <- area_summary$Total_ha / 100

area_summary

# Total valid pixels assessed
total_pixels <- global(!is.na(LandtrendR$yod), "sum")[1,1]

total_pixels


# Percentage of total area experiencing change
area_summary$Percentage_Change <- 
  (area_summary$Change_ha / area_summary$Total_ha) * 100

# Difference in changed area (ha)
difference_ha <- change_area_2 - change_area_1

# Relative increase in changed area (%)
relative_increase <- 
  (difference_ha / change_area_1) * 100

difference_ha
relative_increase
#change_20182022_poly <- as.polygons(change_20182022)
#plot(change_20182022_poly)
#filepath <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/change_20182022_vector.kml"
#writeVector(change_20182022_poly,filepath , overwrite = TRUE)
#filepath <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/change_20182022.tif"
#terra::writeRaster(change_20182022,filepath, overwrite = TRUE)


writeRaster(your_spatraster, "path/to/output.tif", overwrite = TRUE)

#🟧🟨🟧🟨🟧🟨🟧🟨 Local Moran's I values 🟨🟧🟨🟧🟨🟧🟨🟧🟨
set.seed(123)
MI13<- autocor(as.numeric(change_20132017), method = "moran", w=3, global = FALSE)
#filepath <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/MI13.tif"
#terra::writeRaster(MI13,filepath , overwrite = TRUE)

set.seed(123)
MI18<- autocor(as.numeric(change_20182022), method = "moran", w=3, global = FALSE)
#filepath <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/MI18.tif"
#terra::writeRaster(MI18,filepath , overwrite = TRUE)
#MI18 <- rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/MI18.tif")
MI18 <- as.numeric(MI18)
MI18_poly <- as.polygons(MI18)

#MI18 <- rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/MI18.tif")


R.version
citation("terra")
library(terra)
#🟨🟧🟨🟧🟨🟧🟨🟧🟨 Build full datset dataframe 

df <- data.frame(
  Moran_I = c(values(MI13, na.rm = TRUE),
              values(MI18, na.rm = TRUE)),
  Period = factor(
    rep(c("2013–2017", "2018–2022"),
        c(ncell(MI13) - global(is.na(MI13), "sum")[1,1],
          ncell(MI18) - global(is.na(MI18), "sum")[1,1]))
  ), 
  Change = c(values(change_20132017, na.rm = TRUE),
             values(change_20182022, na.rm = TRUE)), 
  Y = c(crds(change_20132017, na.rm = TRUE)[,2],
        crds(change_20182022, na.rm = TRUE)[,2]),
  X = c(crds(change_20132017, na.rm = TRUE)[,1],
        crds(change_20182022, na.rm = TRUE)[,1])
)

#🟧🟨🟧🟨🟧🟨🟧🟨 Clustering over time 

hist(df$Moran_I)
hist(log1p(df$Moran_I))
hist(log(df$Moran_I))
df$Change[df$Change == "FALSE"] <- "Non_Browned"
df$Change[df$Change == "TRUE"] <- "Browned"
df$Change <- as.factor(df$Change)
df$Change  <- factor(df$Change, levels = c("Non_Browned", "Browned"))


set.seed(123)
mod_1 <- lm(log1p(Moran_I)~ Period*Change, data=df)
summary(mod_1)

# Predicted values for each combination
newdat <- expand.grid(
  Period = levels(df$Period),
  Change = levels(df$Change)
)
library(emmeans)
emm <- emmeans(mod_1, ~ Period * Change)

#contrast(emm, interaction = "revpairwise")

contrast(emm, "revpairwise", by = "Change")   # period change within each group
#contrast(emm, "revpairwise", by = "Period")   # browned vs non-browned within each period

ct <- contrast(emm, "revpairwise", by = "Change")
summary(ct, infer = TRUE)   # estimates, CIs, t, p

newdata <- expand.grid(
  Period = levels(df$Period),
  Change = levels(df$Change)
)
newdat$log1p_pred <- predict(mod_1, newdata = newdat)
newdat$Moran_I_pred <- exp(newdat$log1p_pred) - 1
newdat

confint(mod_1)
emm <- emmeans(mod_1, ~ Period * Change)
summary(emm)                       # log1p scale with CIs
cbind(summary(emm)[, 1:2], lapply(summary(emm)[, c("emmean","lower.CL","upper.CL")], expm1))


# ---- 0. Rows dropped by log1p? ----
sum(is.na(log1p(df$Moran_I)))       # NaN where Moran_I <= -1
nobs(mod_1) == nrow(df)             # FALSE means rows were dropped

# ---- 1. Residuals and fitted values ----
d <- model.frame(mod_1)             # only rows actually used
d$res <- resid(mod_1)
d$fit <- fitted(mod_1)

# Subsample for plotting (millions of points is slow and uninformative)
set.seed(123)
ds <- d[sample(nrow(d), min(50000, nrow(d))), ]

dev.cur() 
dev.off()     
graphics.off()
citation("stats")

# ---- 2. Normality of residuals ----
qqnorm(ds$res); qqline(ds$res)
hist(ds$res, breaks = 100)

# ---- 3. Constant variance across groups ----
ggplot(ds, aes(Period:Change, res)) + geom_boxplot()
aggregate(res ~ Period + Change, data = d, FUN = sd)

# ---- 4. Residuals vs fitted ----
plot(ds$fit, ds$res, pch = ".", xlab = "Fitted", ylab = "Residual")
abline(h = 0, col = "red")

# ---- 5. Group sizes (check balance) ----
table(d$Period, d$Change)

p3 <-ggplot(df, aes(x = Change, y = Moran_I, fill = Period)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values=c("#a6bddb", "#1c9099")) + 
  theme_bw() + 
  xlab("Category") +
  ylab("Local Moran's I") +
  theme(text = element_text(size = 20))+
  scale_y_continuous(breaks=seq(-1, 7, 1))

p3
ggsave(
  "p3_local_morans_I.pdf", p3,
  width = 8, height = 6, device = cairo_pdf
)


library(broom)
library(dplyr)

# Coefficient table
coef_tbl <- broom::tidy(mod_1, conf.int = TRUE) |>
  transmute(
    Term       = term,
    Estimate   = round(estimate, 3),
    `Std. Error` = round(std.error, 3),
    `95% CI`   = paste0("[", round(conf.low, 3), ", ", round(conf.high, 3), "]"),
    `t value`  = round(statistic, 2),
    `p value`  = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
  )

# Tidy up term names
coef_tbl$Term <- c(
  "Intercept (2013–2017, reference category)",
  "Period 2018–2022",
  "Change: Browned",
  "Period 2018–2022 × Browned"
)

# Model fit statistics
fit_tbl <- broom::glance(mod_1) |>
  transmute(
    `R²`          = round(r.squared, 3),
    `Adjusted R²` = round(adj.r.squared, 3),
    `Residual SE` = round(sigma, 3),
    `F statistic` = round(statistic, 0),
    `df`          = paste(df, df.residual, sep = ", "),
    `p value`     = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
  )

coef_tbl
fit_tbl

# Export
write.csv(coef_tbl, "mod_1_coefficients.csv", row.names = FALSE)
write.csv(fit_tbl,  "mod_1_fit.csv",          row.names = FALSE)








#🟨🟧🟨🟧🟨🟧🟨🟧🟨 Environmental predictors

file_path_Elev <- "/Users/catherinemercer/Downloads/elevation.tif"
elevation <- terra::rast(file_path_Elev)
crs(elevation)
pts <- terra::vect(
  df,
  geom = c("X", "Y"),
  crs = "EPSG:4326"
)
ext <- terra::extract(elevation, pts)
df$Elevation <- ext[, 2]

#SLOPE ⛰
Slope <- terra::rast("/Users/catherinemercer/Downloads/slope.tif")
pts <- terra::vect(
  df,
  geom = c("X", "Y"),
  crs = "EPSG:4326"
)
ext <- terra::extract(Slope, pts)
df$Slope <- ext[, 2]

#EASTNESS and NORTHNESS ⬅️

Aspect <- terra::rast("/Users/catherinemercer/Downloads/aspect.tif")
aspect_radians <- Aspect * pi / 180
eastness <- sin(aspect_radians)
ext <- terra::extract(eastness, pts)
df$eastness <- ext[, 2]
northness <- cos(aspect_radians)
ext <- terra::extract(northness, pts)
df$northness <- ext[, 2]

#  TWI ️🌧️
TWI <- terra::rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/TWI_renameNA_inR_181125.tif")
plot(TWI)
ext <- terra::extract(TWI, pts)
df$TopWaterIndex <- ext[, 2]

# Mean curvature ️⤵️
Curve_mean <- terra::rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter1/Resubmission/Topographic/Curvature/MeanCurvature_twnclip_280925.tif")
ext <- terra::extract(Curve_mean, pts)
df$Curvature <- ext[, 2]


# 🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Random subset
set.seed(123)
balanced_data <- df %>%
  group_by(Period, Change) %>% 
  slice_sample(n = 10000, replace = FALSE) %>% # Use replace = TRUE if any combination has fewer than 10k rows
  ungroup()


# filepath <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/random_pts_110826.csv"
# write.csv(balanced_data, filepath)


# 🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Export to retrieve climate data 
subset <- read.csv("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Resubmission/Dieback/random_pts_110826.csv")
levels(subset$Period)

set.seed(123)
subset$Period <- as.factor(subset$Period)
mod_1 <- gam(Change ~ s(Y,by=Period) +s(X,by=Period)+ s(Elevation,by=Period) + s(Slope, by=Period) + s(eastness,by=Period) + s(northness,by=Period) + s(TopWaterIndex , by=Period) + s(Curvature , by=Period), data=subset, family=binomial) 

gam.check(mod_1)# Gam check
names(balanced_data)

summary(mod_1)
plot(mod_1, pages=1, shift=coef(mod_1)[1], trans=plogis) 

# Print table 
create_modelSummary(list(mod_1), method_expTransform="simple", dat = balanced_data)


# 🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧   


# 1. Fit the null model (dropping the 'by=Period' conditional smooths)
set.seed(123)
mod_null <- gam(Change ~ s(Y) + s(X)+ s(Elevation) + s(Slope) + s(eastness) + s(northness) + s(TopWaterIndex) + s(Curvature), data=subset, family=binomial) # REMOVE Y????

# 2. Compare using an Analysis of Deviance (Chisq test for binomial models)
anova(mod_null, mod_1, test = "Chisq")

# 🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Plot indiviual plots 
#Back transforming questions? 

Y,by=Period) + s(Elevation,by=Period) + s(Slope, by=Period) + s(eastness,by=Period) + s(northness,by=Period) + s(TopWaterIndex , by=Period) + s(Curvature , by=Period), data=subset, family=binomial) 

sm_Lattitude <- smooth_estimates(mod_1, select = "s(Latt)", partial_match = TRUE)
graphics.off()



sm_Curvature_plot<- sm_Curvature|>
  add_confint() |>
  ggplot(aes(y = .estimate, x = Curvature, colour = Period)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
              alpha = 0.2,
              fill = "grey",
              linetype = 0) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = c(
      "2013–2017" = "#1c9099",
      "2018–2022" = "#f03b20"
    )
  ) +
  labs(
    y = "Partial effect",
    x = "Curvature",
    colour = "Period"
  ) +
  theme_classic() +
  theme(text = element_text(size = 20), legend.position = "none")

sm_Y_plot
sm_X_plot
sm_Elevation_plot
sm_Slope_plot
sm_Eastness_plot
sm_Northness_plot
sm_TopWaterIndex
sm_Curvature
library(patchwork)
sm_Y_plot + sm_X_plot + sm_Elevation_plot +  sm_Slope_plot + sm_Eastness_plot + sm_Northness_plot + sm_TopWaterIndex_plot + sm_Curvature_plot

combined_plot <- (
  sm_Y_plot +
    sm_X_plot +
    sm_Elevation_plot +
    sm_Slope_plot +
    sm_Eastness_plot +
    sm_TopWaterIndex +
    sm_Curvature
) +
  plot_layout(ncol = 4)

combined_plot


# 🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Back transforming estimates and plotting

levels(balanced_data$Period)

newdata <- data.frame(
  X             = seq(min(balanced_data$X), max(balanced_data$X), length.out = 100),
  Period        = "2018–2022", 
  Y             = mean(balanced_data$Y, na.rm = TRUE),
  Elevation     = mean(balanced_data$Elevation, na.rm = TRUE),
  Slope         = mean(balanced_data$Slope, na.rm = TRUE),
  eastness      = mean(balanced_data$eastness, na.rm = TRUE),
  northness     = mean(balanced_data$northness, na.rm = TRUE),
  TopWaterIndex = mean(balanced_data$TopWaterIndex, na.rm = TRUE),
  Curvature     = mean(balanced_data$Curvature, na.rm = TRUE)
)
# Add another coloum of est_predictor_time period

# 1. Predict on the LINK (logit) scale to get standard errors
pred_link <- predict(mod_1, newdata = newdata, type = "link", se.fit = TRUE)

# 2. Calculate 95% Confidence Intervals on the link scale
lower_link <- pred_link$fit - (1.96 * pred_link$se.fit)
upper_link <- pred_link$fit + (1.96 * pred_link$se.fit)

# 3. Back-transform everything to the RESPONSE (probability) scale
X_2018_2022 <- data.frame(
  X            = new_data$X,
  Period       = new_data$Period,
  Probability  = plogis(pred_link$fit),
  Lower_95_CI  = plogis(lower_link),
  Upper_95_CI  = plogis(upper_link),
  Period = "2018–2022", 
  Predictor= "X"
)

head(final_results)


library(ggplot2)

final_results<- cbind(X_2013_2017,)

ggplot(X_2018_2022, aes(x = X, y = Probability)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Lower_95_CI, ymax = Upper_95_CI), alpha = 0.2, fill = "blue") +
  labs(x = "Variable X", y = "Probability of Change", title = "Effect of X during Period 1") +
  theme_minimal() +
  ylim(0, 1) 






# 🟨🟧🟨🟧🟨🟧🟨🟧🟨🟧  Proportion and frequncy distribution 


LandtrendR <- rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/LandtrendR/lt-loss_map_060325.tif")
plot(LandtrendR)
change_20132017 <- terra::rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/LandtrendR/04_03)26/Change_20132017.tif")
change_20132017_df <- as.data.frame(change_20132017)

nochange_20132017 <- terra::rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/LandtrendR/04_03)26/20132017_nochange.tif") # Looks okay
nochange_20132017_df <- as.data.frame(nochange_20132017)

change_20182022 <- terra::rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/LandtrendR/04_03)26/20182022_change.tif")# Looks okay
change_20182022_df <- as.data.frame(change_20182022)

nochange_20182022 <- terra::rast("/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/LandtrendR/04_03)26/20182022_nochange.tif")# Looks okay
nochange_20182022_df<- as.data.frame(nochange_20182022)



change_20132017_df$Type <- "Change"
nochange_20132017_df$Type <- "No Change"
change_20182022_df$Type <- "Change"
nochange_20182022_df$Type <- "No Change"

change_20132017_df$TimePeriod <- "2013-2017"
nochange_20132017_df$TimePeriod <- "2013-2017"
change_20182022_df$TimePeriod <- "2018-2022"
nochange_20182022_df$TimePeriod <- "2018-2022"

year_20132017 <- rbind(change_20132017_df,nochange_20132017_df)
year_combined<- rbind(year_20132017,change_20182022_df)
year_all <- rbind(year_combined,nochange_20182022_df)
unique(year_all)

library(tidyverse)
names(year_all)

Proportion <- year_all %>%
  count(TimePeriod, Type) %>%
  group_by(TimePeriod) %>%
  mutate(freq = n / sum(n))

Proportion_df <- as.data.frame(Proportion)
Proportion_df$Type <- as.factor(Proportion_df$Type)
Proportion_df$Type  <- factor(Proportion_df$Type , levels = c("No Change", "Change"))
names(colors) <- levels(Proportion_df$Type)



p1 <- Proportion_df %>% 
  ggplot(aes(x = TimePeriod, y = freq, fill = Type)) + 
  geom_col(position = "dodge", colour = "black") +
  theme_bw() + 
  xlab("Time Period") +
  ylab("Proportion of total pixels") +
  guides(color = none)+
  scale_fill_manual(
    values =  c("#addd8e", "#fdae6b"),
    labels = c("Non-browned", "Browned")
  )+
  theme(text = element_text(size = 20))




library(tidyverse)

# Counts of browned / total pixels per time period
Proportion_df <- year_all %>%
  count(TimePeriod, Type) %>%
  group_by(TimePeriod) %>%
  mutate(total = sum(n),
         freq  = n / total) %>%
  ungroup() %>%
  filter(Type == "Change") %>%
  rowwise() %>%
  mutate(
    ci_low  = prop.test(n, total)$conf.int[1],
    ci_high = prop.test(n, total)$conf.int[2]
  ) %>%
  ungroup()

p1 <- Proportion_df %>%
  ggplot(aes(x = TimePeriod, y = freq)) +
  geom_col(fill = "#fdae6b", colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8) +
  theme_bw() +
  xlab("Time Period") +
  ylab("Proportion of pixels browned") +
  theme(text = element_text(size = 20))

p1



# Raster properties should match within a period
compareGeom(change_20132017, nochange_20132017)
compareGeom(change_20182022, nochange_20182022)

# What values do the rasters actually hold, and how many NAs?
freq(change_20132017)
freq(nochange_20132017)
freq(change_20182022)
freq(nochange_20182022)

# Total cells vs non-NA cells
ncell(change_20132017)
global(change_20132017, "notNA")
global(nochange_20132017, "notNA")

# Overlap between change and no-change (should be 0)
global(!is.na(change_20132017) & !is.na(nochange_20132017), "sum")
global(!is.na(change_20182022) & !is.na(nochange_20182022), "sum")

# Do the two periods share the same footprint?
global(!is.na(change_20132017) | !is.na(nochange_20132017), "sum")
global(!is.na(change_20182022) | !is.na(nochange_20182022), "sum")


# Calculating pixel area
pixel_area <- 30*30
metres2 <- pixel_area *1379342 #see above
km2 <- metres2/1000000

#### MagnitudeChange 🟧🟧🟧🟧🟨🟨🟨🟨
# Proportion of low, median, high magnitude pixels 
# Only need to do for change pixels - no change should have no magntidue

change_20132017_magnitude <- terra::mask(
  LandtrendR$mag,
  change_20132017
)
 # looks similar to above
change_20132017_magnitude_df <- as.data.frame(change_20132017_magnitude)
change_20132017_magnitude_df$TimePeriod <- "2013-2017"
# %>% rename("2013-2017" = 'mag')

change_20182022_magnitude <- terra::mask(
  LandtrendR$mag,
  change_20182022
)
# looks similar to above
library(tidyverse)
change_20182022_magnitude_df <- as.data.frame(change_20182022_magnitude)
change_20182022_magnitude_df$TimePeriod <- "2018-2022"


# Frequency diagrams  🌸🌸🌸🌸🌸

mag_change_joined <- rbind(change_20182022_magnitude_df,change_20132017_magnitude_df)
summary(mag_change_joined)

mag_change_joined$TimePeriod <- as.factor(mag_change_joined$TimePeriod )


library(tidyverse)
library(hrbrthemes)
library(viridis)

p2 <- ggplot(data=mag_change_joined, aes(x=mag, fill=TimePeriod)) +
  geom_density(adjust=1.5, alpha=.4)+
  scale_fill_manual("Time Period", values = c("2018-2022"= "red", "2013-2017" = "yellow"))+
  theme_bw() + 
  xlab("Magnitude of Browning") +
  ylab("Proportion of Browned Pixels") +
  theme(text = element_text(size = 20))

p1 <- p1 + theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) + theme(text = element_text(size = 20))
p2 <- p2 + theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) + theme(text = element_text(size = 20))
p3 <- p3 + theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) + theme(text = element_text(size = 20))

p1 + p2 + p3


