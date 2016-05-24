## ---- echo=FALSE---------------------------------------------------------
library(clifro)
library(pander)
surfaceWind.dt = new("cfDatatype"
    , dt_name = "Wind"
    , dt_type = "Surface wind"
    , dt_sel_option_names = list("9amWind")
    , dt_sel_combo_name = "knots"
    , dt_param = structure("ls_sfw,1,2,3,4,5", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("132", "knots"), .Names = c("prm4", "prm5")))
    , dt_selected_options = list(c(4, 5))
    , dt_option_length = 5
)

menu.opts = function(title, options){
  cat(paste(title, "",
              paste(seq_along(options), options, sep = ": ", 
                    collapse = "\n"), sep = "\n"))
}

## ---- eval=FALSE---------------------------------------------------------
#  surfaceWind.dt = cf_datatype()
#  
#  # If you prefer pointing and clicking - turn the graphics option on:
#  surfaceWind.dt = cf_datatype(graphics = TRUE)

## ---- echo=FALSE---------------------------------------------------------
menu.opts("Daily and Hourly Observations", 
          c("Combined Observations", "Wind", "Precipitation", 
                           "Temperature and Humidity", "Sunshine and Radiation", 
                           "Weather", "Pressure", "Clouds", 
                           "Evaporation / soil moisture"))

## ---- echo=FALSE---------------------------------------------------------
menu.opts("Wind", c("Surface wind", "Max Gust"))

## ---- echo=FALSE---------------------------------------------------------
menu.opts("Surface wind options", c("WindRun", "HlyWind", "3HlyWind", "9amWind")
          )

## ---- echo=FALSE---------------------------------------------------------
menu.opts("Choose another option?", c("yes", "no"))

## ---- echo=FALSE---------------------------------------------------------
menu.opts("Units", c("m/s", "km/hr", "knots"))

## ------------------------------------------------------------------------
surfaceWind.dt

## ---- eval = FALSE-------------------------------------------------------
#  surfaceWind.dt = cf_datatype(2, 1, 4, 3)
#  surfaceWind.dt

## ---- echo = FALSE-------------------------------------------------------
surfaceWind.dt

## ---- echo = FALSE-------------------------------------------------------
surfaceWind.dt = new("cfDatatype"
    , dt_name = "Wind"
    , dt_type = "Surface wind"
    , dt_sel_option_names = list(c("HlyWind", "9amWind"))
    , dt_sel_combo_name = "knots"
    , dt_param = structure("ls_sfw,1,2,3,4,5", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("134", "132", "knots"), .Names = c("prm2", "prm4", 
"prm5")))
    , dt_selected_options = list(c(2, 4, 5))
    , dt_option_length = 5
)

rainfall.dt = new("cfDatatype"
    , dt_name = "Precipitation"
    , dt_type = "Rain (fixed periods)"
    , dt_sel_option_names = list(c("Daily ", "Hourly"))
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_ra,1,2,3,4", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("181", "182"), .Names = c("prm1", "prm2")))
    , dt_selected_options = list(c(1, 2))
    , dt_option_length = 4
)

lightning.dt = new("cfDatatype"
    , dt_name = "Weather"
    , dt_type = "Lightning"
    , dt_sel_option_names = list("Ltng")
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_light,1", .Names = "dt1")
    , dt_sel_option_params = list(structure("271", .Names = "prm1"))
    , dt_selected_options = list(1)
    , dt_option_length = 1
)

temperatureExtremes.dt = new("cfDatatype"
    , dt_name = "Temperature and Humidity"
    , dt_type = "Max_min_temp"
    , dt_sel_option_names = list(c("DlyGrass", "HlyGrass"))
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_mxmn,1,2,3,4,5,6", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("202", "204"), .Names = c("prm5", "prm6")))
    , dt_selected_options = list(c(5, 6))
    , dt_option_length = 6
)

## ---- eval = FALSE-------------------------------------------------------
#  surfaceWind.dt = cf_datatype(2, 1, c(2, 4), 3)
#  surfaceWind.dt

## ---- echo = FALSE-------------------------------------------------------
surfaceWind.dt

## ---- eval=FALSE---------------------------------------------------------
#  # Hourly and 9am surface wind (knots)
#  surfaceWind.dt = cf_datatype(2, 1, c(2, 4), 3)
#  surfaceWind.dt

## ---- echo = FALSE-------------------------------------------------------
surfaceWind.dt

## ---- eval = FALSE-------------------------------------------------------
#  # Hourly and daily rainfall
#  rainfall.dt = cf_datatype(3, 1, c(1, 2))
#  rainfall.dt

## ---- echo = FALSE-------------------------------------------------------
rainfall.dt

## ---- eval = FALSE-------------------------------------------------------
#  # Hourly counts of lightning flashes
#  lightning.dt = cf_datatype(6, 1, 1)
#  lightning.dt

## ---- echo = FALSE-------------------------------------------------------
lightning.dt

## ---- eval = FALSE-------------------------------------------------------
#  # Daily and hourly grass temperature extremes
#  temperatureExtremes.dt = cf_datatype(4, 2, c(5, 6))
#  temperatureExtremes.dt
#  
#  # Note: only the surface wind datatype requires combo options

## ---- echo = FALSE-------------------------------------------------------
temperatureExtremes.dt

## ---- echo = FALSE, results = "asis"-------------------------------------
d = data.frame(Menu = c("First selection", "Second selection", 
                        "Third selection(s)", "combo box options"),
               `Surface wind` = c(2, 1, "2 & 4", 3),
               Rainfall = c(3, 1, "1 & 2", NA),
               Lightning = c(6, 1, 1, NA),
               Temperature = c(4, 2, "5 & 6", NA))
pandoc.table(d, style = "simple")

## ---- echo = FALSE-------------------------------------------------------
query1.dt = new("cfDatatype"
    , dt_name = c("Wind", "Precipitation", "Weather", "Temperature and Humidity"
)
    , dt_type = c("Surface wind", "Rain (fixed periods)", "Lightning", "Max_min_temp"
)
    , dt_sel_option_names = list(c("HlyWind", "9amWind"), c("Daily ", "Hourly"), "Ltng", 
    c("DlyGrass", "HlyGrass"))
    , dt_sel_combo_name = c("knots", NA, NA, NA)
    , dt_param = structure(c("ls_sfw,1,2,3,4,5", "ls_ra,6,7,8,9", "ls_light,10", 
"ls_mxmn,11,12,13,14,15,16"), .Names = c("dt1", "dt2", "dt3", 
"dt4"))
    , dt_sel_option_params = list(structure(c("134", "132", "knots"), .Names = c("prm2", "prm4", 
"prm5")), structure(c("181", "182"), .Names = c("prm6", "prm7"
)), structure("271", .Names = "prm10"), structure(c("202", "204"
), .Names = c("prm15", "prm16")))
    , dt_selected_options = list(c(2, 4, 5), c(1, 2), 1, c(5, 6))
    , dt_option_length = c(5, 4, 1, 6)
)

## ---- tidy = FALSE, eval = FALSE-----------------------------------------
#  query1.dt = cf_datatype(c(2, 3, 6, 4),
#                          c(1, 1, 1, 2),
#                          list(c(2, 4), c(1, 2), 1, c(5, 6)),
#                          c(3, NA, NA, NA))
#  query1.dt

## ---- echo = FALSE-------------------------------------------------------
query1.dt

## ------------------------------------------------------------------------
query1.dt = surfaceWind.dt + rainfall.dt + lightning.dt + 
  temperatureExtremes.dt
query1.dt

## ---- eval=FALSE---------------------------------------------------------
#  # To add another datatype using the menu:
#  query1.dt + cf_datatype()
#  
#  # Is equivalent to:
#  query1.dt + cf_datatype(NA, NA, NA, NA)
#  
#  # Therefore is equivalent to adding a column of NA's to the above table:
#  query1.dt = cf_datatype(c(2, 3, 6, 4, NA),
#                                c(1, 1, 1, 2, NA),
#                                list(c(2, 4), c(1, 2), 1, c(5, 6), NA),
#                                c(3, NA, NA, NA, NA))
#  
#  # Half an unknown wind datatype i.e. we know first selection = 2 but nothing
#  # further:
#  rain.dt = cf_datatype(2) # Or cf_datatype(2, NA, NA, NA)

