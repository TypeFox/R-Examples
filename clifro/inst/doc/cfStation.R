## ---- echo = FALSE-------------------------------------------------------
library(clifro)

## ---- eval = FALSE-------------------------------------------------------
#  lake.tekapo.st = cf_station(12709, 35567, 39557, 4630, 24945, 4616, 4602)
#  lake.tekapo.st[, c("name", "agent", "start", "end", "open")]

## ---- eval = FALSE-------------------------------------------------------
#  added.stations.st = lake.tekapo.st +
#    cf_station() +
#    cf_find_station("lighthouse", status = "all")
#  added.stations.st[, c("name", "agent", "start", "end", "open")]

## ---- echo = FALSE-------------------------------------------------------
auckland.df = dget(system.file("extdata", "auckStations", package = "clifro"))
auckland.st = new("cfStation", auckland.df)

## ---- eval = FALSE-------------------------------------------------------
#  # Conduct the search
#  auckland.st = cf_find_station("auckland", search = "region", status = "all")

## ---- message = FALSE----------------------------------------------------
library(ggmap)

# Add a column to colour the open and closed stations
auckland.st$colour = factor(auckland.st$open, labels = c("Closed", "Open"))

# Coerce to a data.frame and reverse the rows so the open stations get plotted 
# on top of the closed stations
auckland.df = as(auckland.st, "data.frame")[nrow(auckland.st):1, ]

# Obtain the map of the greater Auckland suitably scaled to fit the stations
auckland.map = ggmap(get_map("Auckland", maptype = "hybrid", zoom = 8))

# Plot the resulting map with the stations and station density
auckland.map %+% auckland.df + 
  stat_density2d(aes(colour = colour), alpha = .8) +
  geom_point(aes(colour = colour), alpha = .5) +
  scale_colour_discrete("Status", c("Closed", "Open")) +
  theme(legend.title = element_text(face = "bold"))

