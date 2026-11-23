# =============================================================================
# Taiwan alpine greening: LandTrendr NDVI gain (TWD97 / TM2 zone 121, EPSG:3826)
# =============================================================================
library(terra)
library(tidyterra)
library(dplyr)
library(ggplot2)
library(patchwork)
library(mgcv)
library(gratia)

# ---- 0. Load data ------------------------------------------------------------
base <- "/Volumes/Mercer_PhD_2/PhD_Data/Chapter2/Revised_chpt3_alpine/"
Shrub_Herb <- terra::vect(paste0(base, "Shrub_and_herbaceous.shp"))
r <- rast(paste0(base, "NDVI_gain_1990_2025.tif"))

nlyr(r)                       # must be 5
names(r) <- c("yod", "mag", "elev", "slope", "aspect")
summary(r)                    # yod: 0 or a year; mag: 0 or positive
crs(r, describe = TRUE)       # expect EPSG:3826 (metres)

e <- as.data.frame(r)
summary(e)

# ---- 1. Data frame -----------------------------------------------------------
# Study period: only gains detected 2000-2025 are analysed.
# Pixels with a gain detected before 2000 are dropped (not coded as "not greened").
yr_min <- 2000
yr_max <- 2025

df <- as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
  filter(slope > 3) |>                       # aspect meaningless on flat ground
  filter(yod == 0 | (yod >= yr_min & yod <= yr_max)) |>
  mutate(
    greened = yod > 0,                       # TRUE if a gain segment was detected
    north   = cos(aspect * pi / 180),        # +1 = north-facing
    east    = sin(aspect * pi / 180)         # +1 = east-facing
  )

names(df)
table(df$greened)
range(df$yod[df$yod > 0])                    # effective detection window

yearly <- df |>
  filter(greened) |>
  group_by(yod) |>
  summarise(n_pixels = n(),
            area_km2 = n() * 900 / 1e6,
            mean_mag = mean(mag),
            sd_mag   = sd(mag),
            .groups  = "drop") |>
  mutate(pct_of_area = 100 * n_pixels / nrow(df))
print(yearly, n = Inf)

# ---- 2. GAMs -----------------------------------------------------------------
set.seed(123)
sub <- df |> slice_sample(n = 20000)

hist(df$north); hist(df$east)

# Presence of greening
set.seed(123)
m1 <- bam(greened ~ s(elev) + s(slope) + s(north) + s(east) + s(x),
          family = binomial, data = sub)

# Magnitude where greening occurred
set.seed(123)
m2 <- bam(mag ~ s(elev) + s(slope) + s(north) + s(east) + s(x),
          data = filter(sub, greened))
summary(m1); summary(m2)

# Quick looks
plot(m1, pages = 1, pch = 16, cex = 0.3)
plot(m2, pages = 1, pch = 16, cex = 0.3)
draw(m1); draw(m2)

# ---- 3. Shared plot style and size (matches the m2_x_plot style) ----------------
theme_ms <- function() {
  theme_classic() +
    theme(text = element_text(size = 22), legend.position = "none")
}
plot_w <- 7   # inches, used for every saved plot
plot_h <- 6

# ---- 3b. Individual smooth plots (helper) ---------------------------------------
smooth_plot <- function(model, term, xlab, ylim = NULL, xscale = 1) {
  est <- smooth_estimates(model, select = paste0("s(", term, ")"),
                          partial_match = TRUE) |>
    add_confint()
  p <- ggplot(est, aes(x = .data[[term]] / xscale, y = .estimate)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
                alpha = 0.2, fill = "grey", linetype = 0) +
    geom_line(linewidth = 1) +
    labs(x = xlab, y = "Partial effect") +
    theme_ms()
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# m1: occurrence (log-odds)
set.seed(123)
m1_elev_plot  <- smooth_plot(m1, "elev",  "Elevation (m)",        c(-4, 4))
set.seed(123)
m1_slope_plot <- smooth_plot(m1, "slope", "Slope (°)",            c(-4, 4))
set.seed(123)
m1_north_plot <- smooth_plot(m1, "north", "Northness",            c(-4, 4))
set.seed(123)
m1_east_plot  <- smooth_plot(m1, "east",  "Eastness",             c(-4, 4))
set.seed(123)
m1_x_plot     <- smooth_plot(m1, "x",     "Easting (km, TWD97)",  c(-4, 4), xscale = 1000)

m1_all <- (m1_elev_plot + m1_slope_plot + m1_north_plot) /
  (m1_east_plot + m1_x_plot + plot_spacer())
m1_all

# m2: magnitude
set.seed(123)
m2_elev_plot  <- smooth_plot(m2, "elev",  "Elevation (m)",        c(-40, 40))
set.seed(123)
m2_slope_plot <- smooth_plot(m2, "slope", "Slope (°)",            c(-40, 40))
set.seed(123)
m2_north_plot <- smooth_plot(m2, "north", "Northness",            c(-40, 40))
set.seed(123)
m2_east_plot  <- smooth_plot(m2, "east",  "Eastness",             c(-40, 40))
set.seed(123)
m2_x_plot     <- smooth_plot(m2, "x",     "Easting (km, TWD97)",  c(-40, 40), xscale = 1000)

m2_all <- (m2_elev_plot + m2_slope_plot + m2_north_plot) /
  (m2_east_plot + m2_x_plot + plot_spacer())
m2_all

# ---- 4. LandTrendr timing plots ----------------------------------------------
# Be cautious near the ends of the record and at sensor transitions (L7, L8).
# Reusable x-axis settings: every 2 years, labels tilted 45 degrees
every_2_years <- list(
  scale_x_continuous(breaks = seq(yr_min, yr_max, by = 2),
                     limits = c(yr_min - 0.5, yr_max + 0.5)),
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

x <- ggplot(yearly, aes(yod, area_km2)) +
  geom_col(fill = "olivedrab", alpha = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
              se = TRUE, level = 0.95,
              colour = "black", fill = "grey40", alpha = 0.25,
              linewidth = 0.8) +
  labs(x = "Year of detection", y = expression("Area greened (km"^2*")")) +
  theme_ms() +
  every_2_years

p2 <- yearly |> mutate(cum_km2 = cumsum(area_km2)) |>
  ggplot(aes(yod, cum_km2)) + geom_line() + geom_point() +
  labs(x = "Year", y = expression("Cumulative greening area (km"^2*")")) +
  theme_ms() +
  every_2_years

p1 + p2

ggplot(yearly, aes(yod, mean_mag)) +
  geom_point() + geom_line() +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick",
              fill = "firebrick", alpha = .15) +
  labs(x = "Year of detection", y = "Mean magnitude (NDVI x1000)") +
  theme_ms()
