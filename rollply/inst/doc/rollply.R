## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)

## ---- output = "none", messages = FALSE----------------------------------
library(proj4) 
library(ggplot2) 
library(rgdal)
library(tidyr) 
library(rollply) 
library(plyr) 

## ---- fig.width=7, fig.height=4------------------------------------------
# Download and format data
url <- "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt"
hawaii <- read.table(url)[ ,c(3,4)]
names(hawaii) <- c('date','CO2')
hawaii[hawaii$CO2 < 0, "CO2"] <- NA # mark NAs as such

# Display original trend
CO2.plot <- ggplot(hawaii) + geom_line(aes(date, CO2)) + ylab("CO2 (ppm)")
print(CO2.plot)

## ---- fig.width=7, fig.height=4------------------------------------------
# with smoothed trend
hawaii.smoothed <- rollply(hawaii, ~ date, wdw.size = 1,
                           summarize, CO2.mean = mean(CO2, na.rm = TRUE), )
CO2.plot + geom_line(aes(date, CO2.mean), data = hawaii.smoothed, color = 'red')

## ---- warning=FALSE, results="hide"--------------------------------------
# Download and prepare dataset
# Source and decription: 
# https://publicdata.eu/dataset/listes-des-communes-par-rgions-dpartements-circonscriptions

tmpfile <- tempfile()
url <- paste0('http://www.nosdonnees.fr/wiki/images/b/b5/',
         'EUCircos_Regions_departements_circonscriptions_communes_gps.csv.gz')
download.file(url, destfile = tmpfile)
dat <- read.csv2(tmpfile, stringsAsFactors = FALSE)
file.remove(tmpfile)
dat <- dat[with(dat, latitude != "" | ! grepl(",", latitude)), 
           c('nom_commune', 'latitude', 'longitude')]
colnames(dat) <- c('name', 'latitude', 'longitude')

dat[ ,'name']      <- as.factor(tolower(dat[ ,'name']))
dat[ ,'latitude']  <- as.numeric(dat[ ,'latitude'])
dat[ ,'longitude'] <- as.numeric(dat[ ,'longitude'])

# We use an equirectangular projection to work on true distances
dat <- na.omit(dat)
dat <- data.frame(dat, proj4::project(dat[ ,c('longitude','latitude')],
                                      '+proj=eqc'))
dat <- dat[ ,c('name','x','y')]


## ---- fig.width=7, fig.height=6------------------------------------------
# Visualise distribution of towns
str(dat)
ggplot(dat) + geom_point(aes(x, y), alpha=.1)

## ---- fig.width=7, fig.height=6------------------------------------------

# This is our custom function : it accepts a data frame and a regular 
# expression and returns the number of matches in the column "name", formatted
# withing a data.frame (this plays well with ddply).
how_many_with_name <- function(df, regexp) { 
  data.frame(ker = sum(grepl(regexp, df[ ,'name'])))
}

dat_roll <- rollply(dat, ~ x + y, wdw.size = 1e4, grid_npts = 10e3,
                    how_many_with_name, regexp = "^ker")

ggplot(dat_roll) +
  geom_raster(aes(x, y, fill = ker)) +
  scale_fill_distiller(palette = 'Greys')

## ------------------------------------------------------------------------
dat_roll <- rollply(dat, ~ x + y, wdw.size = 1e4,
                    grid_npts = 10e3, # number of grid points 
                    grid_type = "ahull_fill", # grid type: fills an alpha-hull with the given number of points
                    grid_opts = list(alpha = .05, # shape parameter of the hull
                                     verbose = TRUE), 
                    how_many_with_name, regexp = "^ker")

## ---- fig.width=7, fig.height=6------------------------------------------
ggplot(dat_roll) +
  geom_raster(aes(x, y, fill = ker)) +
  scale_fill_distiller(palette = 'Greys')

## ---- echo=FALSE, fig.width=4*2.5, fig.height=4--------------------------
data(meadow)

# We project lon/lat to UTM zone 11 (southern california)
meadow[ ,c('x','y')] <- proj4::project(as.matrix(meadow[ ,c('lon','lat')]),
                                       "+proj=utm +zone=11 +ellps=WGS84")

ggplot(meadow, aes(x,y)) +
  geom_point(shape='+') +
  ggtitle('Sample points') +
  xlab('UTM X') +
  ylab('UTM Y')

## ------------------------------------------------------------------------
# We request a grid with approximately this number of points:
npts <- 500
base.plot <- ggplot(NULL, aes(x,y)) +
               geom_point(data = meadow, shape='+') +
               xlab('UTM X') +
               ylab('UTM Y') 

grids <- list(identical  = build_grid_identical(meadow[ ,c('x','y')], npts),
              squaretile = build_grid_squaretile(meadow[ ,c('x','y')], npts),
              ahull_crop = build_grid_ahull_crop(meadow[ ,c('x','y')], npts),
              ahull_fill = build_grid_ahull_fill(meadow[ ,c('x','y')], npts))

plot_grid <- function(grid_type) {
  base.plot +
    geom_point(data=grids[[grid_type]]) +
    annotate('text', x = min(meadow$x), y = min(meadow$y),
            label = paste(nrow(grids[[grid_type]]), "points"),
            hjust = 0, vjust = 1)
}


## ---- fig.width=4*2.5, fig.height=4--------------------------------------
plot_grid('identical')

## ---- fig.width=4*2.5, fig.height=4*2.5/4.92-----------------------------
plot_grid('squaretile')

## ---- fig.width=4*2.5, fig.height=4*2.5/4.92-----------------------------
plot_grid('ahull_crop')

## ---- fig.width=4*2.5, fig.height=4*2.5/4.92-----------------------------
plot_grid('ahull_fill')

