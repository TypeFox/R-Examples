## ---- echo=FALSE---------------------------------------------------------
library(clifro)

## ---- eval = FALSE-------------------------------------------------------
#  # Equivalent to searching for status = "open" on CliFro
#  # Note the search string is not case sensitive
#  cf_find_station("takaka", status = "all")

## ---- eval = FALSE-------------------------------------------------------
#  cf_find_station("takaka", status = "open")

## ---- eval = FALSE-------------------------------------------------------
#  cf_find_station("f028", search = "network", status = "all")

## ---- echo = FALSE-------------------------------------------------------
open.queenstown.stations.df = dget(system.file("extdata", "queenStations", package = "clifro"))
open.queenstown.stations = new("cfStation", open.queenstown.stations.df)

## ---- eval = FALSE-------------------------------------------------------
#  # Partial match for the Queenstown region
#  open.queenstown.stations = cf_find_station("queen", search = "region")

## ---- echo = FALSE-------------------------------------------------------
takaka.town.df = structure(list(name = structure(c(2L, 3L, 1L), .Label = c("Takaka Ews", 
"Takaka, Kotinga Road", "Takaka Pohara"), class = "factor"), 
    network = structure(1:3, .Label = c("F02882", "F02884", "F02885"
    ), class = "factor"), agent = c(3788, 3790, 23849), start = structure(c(18273600, 
    520516800, 1020081600), class = c("POSIXct", "POSIXt"), tzone = "Pacific/Auckland"), 
    end = structure(c(1425121200, 1425121200, 1425121200), class = c("POSIXct", 
    "POSIXt"), tzone = "Pacific/Auckland"), open = c(TRUE, TRUE, 
    TRUE), distance = c(2.6, 5.7, 1.6), lat = c(-40.872, -40.845, 
    -40.86364), lon = c(172.809, 172.867, 172.80568)), .Names = c("name", 
"network", "agent", "start", "end", "open", "distance", "lat", 
"lon"), row.names = c(NA, -3L), class = "data.frame")
takaka.town.st = new("cfStation", takaka.town.df)

## ---- eval = FALSE-------------------------------------------------------
#  takaka.town.st = cf_find_station(lat = -40.85, long = 172.8, rad = 10, search = "latlong")
#  takaka.town.st[, -c(8, 9)]

## ---- echo = -1----------------------------------------------------------
takaka.town.st[, -c(8, 9)]

# We may rather order the stations by distance from the township
takaka.town.st[order(takaka.town.st$distance), -c(8, 9)]

## ---- echo = FALSE-------------------------------------------------------
hourly.rain.dt = new("cfDatatype"
    , dt_name = "Precipitation"
    , dt_type = "Rain (fixed periods)"
    , dt_sel_option_names = list("Hourly")
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_ra,1,2,3,4", .Names = "dt1")
    , dt_sel_option_params = list(structure("182", .Names = "prm2"))
    , dt_selected_options = list(2)
    , dt_option_length = 4
)

## ---- eval = FALSE-------------------------------------------------------
#  # Create a clifro datatype for hourly rain
#  hourly.rain.dt = cf_datatype(3, 1, 2)
#  hourly.rain.dt

## ---- echo = FALSE-------------------------------------------------------
hourly.rain.dt

## ---- eval = FALSE-------------------------------------------------------
#  # Conduct the search
#  cf_find_station("takaka", datatype = hourly.rain.dt)

## ---- echo = FALSE-------------------------------------------------------
kaitaia.df = structure(list(name = structure(c(4L, 9L, 3L, 8L, 1L, 6L, 5L, 
7L, 2L), .Label = c("Cape Reinga Aws", "Dargaville 2 Ews", "Kaikohe Aws", 
"Kaitaia Aero Ews", "Kaitaia Ews", "Kerikeri Aerodrome Aws", 
"Kerikeri Ews", "Purerua Aws", "Trounson Cws"), class = "factor"), 
    network = structure(c(2L, 7L, 6L, 9L, 1L, 5L, 3L, 4L, 8L), .Label = c("A42462", 
    "A53026", "A53127", "A53191", "A53295", "A53487", "A53762", 
    "A53987", "A54101"), class = "factor"), agent = c(18183, 
    37131, 1134, 1196, 1002, 37258, 17067, 1056, 25119), start = structure(c(960984000, 
    1244030400, 500727600, 788871600, 788871600, 1214395200, 
    913806000, 1025179200, 1067425200), class = c("POSIXct", 
    "POSIXt"), tzone = "Pacific/Auckland"), end = structure(c(1425294000, 
    1425294000, 1425207600, 1425207600, 1425207600, 1425207600, 
    1424775600, 1423825200, 1423738800), class = c("POSIXct", 
    "POSIXt"), tzone = "Pacific/Auckland"), open = c(TRUE, TRUE, 
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), distance = c(0, 
    0, 0, 0, 0, 0, 0, 0, 0), lat = c(-35.0677, -35.72035, -35.424, 
    -35.129, -34.432, -35.262, -35.135, -35.183, -35.93145), 
    lon = c(173.2874, 173.65153, 173.822, 174.015, 172.682, 173.911, 
    173.262, 173.926, 173.85317)), .Names = c("name", "network", 
"agent", "start", "end", "open", "distance", "lat", "lon"), row.names = c(NA, 
-9L), class = "data.frame")
kaitaia.st = new("cfStation", kaitaia.df)
my.composite.search = takaka.town.st + kaitaia.st

## ---- eval = FALSE-------------------------------------------------------
#  my.composite.search = takaka.town.st + cf_find_station("kaitaia",
#                                                         search = "region",
#                                                         datatype = hourly.rain.dt)
#  my.composite.search

## ---- echo = -1----------------------------------------------------------
my.composite.search

# How long have these stations been open for?
transform(my.composite.search, ndays = round(end - start))[, c(1, 10)]

## ---- echo = FALSE-------------------------------------------------------
all.auckland.df = dget(system.file("extdata", "auckStations", package = "clifro"))
all.auckland.st = new("cfStation", all.auckland.df)

## ----eval = FALSE--------------------------------------------------------
#  # First, search for the stations
#  all.auckland.st = cf_find_station("auckland", search = "region", status = "all")

## ----eval=FALSE----------------------------------------------------------
#  # Then save these as a KML
#  cf_save_kml(all.auckland.st, file_name = "all_auckland_stations")

