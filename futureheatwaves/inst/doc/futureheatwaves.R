## ----echo = FALSE, message = FALSE, warning = FALSE----------------------
library(futureheatwaves)
library(ggplot2)

## ------------------------------------------------------------------------
system.file("extdata/cities.csv", package = "futureheatwaves")

## ----eval = FALSE--------------------------------------------------------
#  # Identify location of example files
#  projection_dir_location <- system.file("extdata/cmip5",
#                                         package = "futureheatwaves")
#  city_file_location <- system.file("extdata/cities.csv",
#                                    package = "futureheatwaves")
#  
#  # Process example files
#  gen_hw_set(out = "example_results",
#             dataFolder = projection_dir_location ,
#             dataDirectories = list("historical" = c(1990, 1999),
#                                          "rcp85" = c(2060, 2079)),
#             citycsv = city_file_location,
#             coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#             tasFilenames = "tas_NorthAmerica_12mo.csv",
#             timeFilenames = "time_NorthAmerica_12mo.csv")

## ------------------------------------------------------------------------
data(hw_datafr)
hw_datafr[1:3, c("hw.number", "mean.temp", "length", "start.date",
                 "mean.temp.quantile", "city")]

## ------------------------------------------------------------------------
system.file("extdata/cmip5", package = "futureheatwaves")

## ------------------------------------------------------------------------
average_mean_temp <- function(hw_datafr){
        out <- mean(hw_datafr$mean.temp)
        return(out)
        }

## ------------------------------------------------------------------------
out <- system.file("extdata/example_results", package = "futureheatwaves")
apply_all_models(out = out, FUN = average_mean_temp)

## ------------------------------------------------------------------------
apply_all_models(out = out, FUN = average_mean_temp, city_specific = TRUE)

## ------------------------------------------------------------------------
excess_deaths <- function(hw_datafr, base_mortality, RR = 1.032){
        hw_datafr <- dplyr::left_join(hw_datafr, base_mortality,
                                      by = "city") %>%
                dplyr::mutate(excess_deaths = base_mort * length * RR)
        out <- sum(hw_datafr$excess_deaths)
        return(out)
}

## ----eval = FALSE--------------------------------------------------------
#  apply_all_models(out = out, FUN = excess_deaths, base_mortality = base_mort)

## ------------------------------------------------------------------------
average_mean_temp

## ----eval = FALSE--------------------------------------------------------
#  data(hw_datafr)

## ----eval = FALSE--------------------------------------------------------
#  gen_hw_set(out = "example_results",
#             dataFolder = projection_dir_location ,
#             dataDirectories = list("historical" = c(1990, 1999),
#                                          "rcp85" = c(2060, 2079)),
#             citycsv = city_file_location,
#             coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#             tasFilenames = "tas_NorthAmerica_12mo.csv",
#             timeFilenames = "time_NorthAmerica_12mo.csv",
#             probThreshold = 0.99)

## ----eval = FALSE--------------------------------------------------------
#  gen_hw_set(out = "example_results",
#             dataFolder = projection_dir_location ,
#             dataDirectories = list("historical" = c(1990, 1999),
#                                          "rcp85" = c(2060, 2079)),
#             citycsv = city_file_location,
#             coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#             tasFilenames = "tas_NorthAmerica_12mo.csv",
#             timeFilenames = "time_NorthAmerica_12mo.csv",
#             thresholdBoundaries = c(2070, 2079))

## ------------------------------------------------------------------------
data(datafr)

## ------------------------------------------------------------------------
head(datafr, 3)
id_of_hws <- IDHeatwavesR(datafr = datafr, threshold = 95, numDays = 2)
head(id_of_hws, 3)

## ---- eval = FALSE-------------------------------------------------------
#  gen_hw_set(out = "example_results",
#             dataFolder = projection_dir_location ,
#             dataDirectories = list("historical" = c(1990, 1999),
#                                         "rcp85" = c(2060, 2079)),
#             citycsv = city_file_location,
#             coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#             tasFilenames = "tas_NorthAmerica_12mo.csv",
#             timeFilenames = "time_NorthAmerica_12mo.csv",
#             IDheatwavesFunction = "IDHeatwavesAlternative")

## ---- eval = FALSE-------------------------------------------------------
#  gen_hw_set(out = "example_results",
#             dataFolder = projection_dir_location ,
#             dataDirectories = list("historical" = c(1990, 1999),
#                                          "rcp85" = c(2060, 2079)),
#             citycsv = city_file_location,
#             coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#             tasFilenames = "tas_NorthAmerica_12mo.csv",
#             timeFilenames = "time_NorthAmerica_12mo.csv",
#             projectionBoundaries = c(2060, 2079))

## ---- eval = FALSE-------------------------------------------------------
#  gen_hw_set(out = "example_results",
#             dataFolder = projection_dir_location ,
#             dataDirectories = list("historical" = c(1990, 1999),
#                                          "rcp85" = c(2060, 2079)),
#             citycsv = city_file_location,
#             coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#             tasFilenames = "tas_NorthAmerica_12mo.csv",
#             timeFilenames = "time_NorthAmerica_12mo.csv",
#             referenceBoundaries = c(1990, 1999))

## ----fig.width = 5, message = FALSE, warning = FALSE---------------------
out <- system.file("extdata/example_results", package = "futureheatwaves")
map_grid(plot_model = "bcc1", out = out)

## ----fig.width = 5, message = FALSE, warning = FALSE, eval = FALSE-------
#  a <- map_grid(plot_model = "bcc1", out = out)
#  a + ggtitle("BCC1 CMIP5 model") + theme_dark()

## ----fig.width = 5, fig.height = 7, eval = FALSE-------------------------
#  library(gridExtra)
#  a <- map_grid(plot_model = "bcc1", out = out)
#  b <- map_grid(plot_model = "ccsm", out = out)
#  grid.arrange(a, b, ncol = 1)

