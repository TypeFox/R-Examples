## ----setup, include=FALSE------------------------------------------------
# Global options
library(knitr)
opts_chunk$set(fig.path="fig/")

## ----install, eval=FALSE-------------------------------------------------
#  install.packages("eurostat")

## ----install2, eval=FALSE------------------------------------------------
#  library(devtools)
#  install_github("ropengov/eurostat")

## ---- echo=TRUE----------------------------------------------------------
library(eurostat)
eurostat.functions <- sort(ls("package:eurostat"))
kable(as.data.frame(eurostat.functions))

## ----get_eurostat_toc, warning=FALSE, message=FALSE----------------------
# Load the package
library(eurostat)
library(rvest)

# Get Eurostat data listing
toc <- get_eurostat_toc()

# Check the first items
library(knitr)
kable(head(toc))

## ----search_eurostat, warning=FALSE, message=FALSE-----------------------
# info about passengers
kable(head(search_eurostat("passenger transport")))

## ----get_id, warning=FALSE, message=FALSE, results='asis'----------------
id <- search_eurostat("Modal split of passenger transport", 
        	             type = "table")$code[1]
print(id)

## ----get_eurostat, warning=FALSE, message=FALSE, results='asis'----------
dat <- get_eurostat(id, time_format = "num")

## ----str_dat, warning=FALSE, message=FALSE-------------------------------
str(dat)

## ----head_dat, warning=FALSE, message=FALSE, results='asis'--------------
kable(head(dat))

## ----get_eurostat_json, warning=FALSE, message=FALSE, results='asis'-----

dat2 <- get_eurostat(id, filters = list(geo = c("EU28", "FI"), lastTimePeriod=1), time_format = "num")
kable(dat2)

## ----json_labels, warning=FALSE, message=FALSE, results='asis'-----------
datl2 <- get_eurostat(id, filters = list(geo = c("EU28", "FI"), 
                                         lastTimePeriod = 1), 
                      type = "label", time_format = "num")
kable(head(datl2))

## ----labels, warning=FALSE, message=FALSE, results='asis'----------------
datl <- label_eurostat(dat)
kable(head(datl))

## ----name_labels---------------------------------------------------------

label_eurostat_vars(names(datl))

## ----vehicle_levels------------------------------------------------------
levels(datl$vehicle)

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
data(efta_countries)
kable(efta_countries)

## ----eu_12---------------------------------------------------------------
dat_eu12 <- subset(datl, geo == "European Union (28 countries)" & time == 2012)
kable(dat_eu12, row.names = FALSE)

## ----eu_vehicles_table---------------------------------------------------
library("tidyr")
dat_eu_0012 <- subset(dat, geo == "EU28" & time %in% 2000:2012)
dat_eu_0012_wide <- spread(dat_eu_0012, vehicle, values)
kable(subset(dat_eu_0012_wide, select = -geo), row.names = FALSE)

## ----trains_table--------------------------------------------------------
dat_trains <- subset(datl, geo %in% c("Austria", "Belgium", "Finland", "Sweden")
                     & time %in% 2000:2012 
                     & vehicle == "Trains")

dat_trains_wide <- spread(dat_trains, geo, values) 
kable(subset(dat_trains_wide, select = -vehicle), row.names = FALSE)

## ----trains_plot, fig.width=8, fig.height=3------------------------------
library(ggplot2)
p <- ggplot(dat_trains, aes(x = time, y = values, colour = geo)) 
p <- p + geom_line()
print(p)

## ----plotGallery, warning=FALSE, message=FALSE---------------------------
library(tidyr)

transports <- spread(subset(dat, time == 2012, select = c(geo, vehicle, values)), vehicle, values)

transports <- na.omit(transports)

# triangle plot
library(plotrix)
triax.plot(transports[, -1], show.grid = TRUE, 
           label.points = TRUE, point.labels = transports$geo, 
           pch = 19)

## ----citation, message=FALSE, eval=TRUE----------------------------------
citation("eurostat")

## ----sessioninfo, message=FALSE, warning=FALSE---------------------------
sessionInfo()

