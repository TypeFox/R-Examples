## ------------------------------------------------------------------------
library(trip)
d <- data.frame(x=1:10, y=rnorm(10), tms=Sys.time() + 1:10, id=gl(2, 5))
coordinates(d) <- ~x+y
## a projection should always be set, is it WGS84 or NAD83 . . .
proj4string(d) <- CRS("+proj=laea +ellps=sphere")
tr <- trip(d, c("tms", "id"))
summary(tr)

## ------------------------------------------------------------------------
plot(tr)
lines(tr)

## ------------------------------------------------------------------------
tg <- rasterize(tr)
plot(tg, col = c("transparent", heat.colors(25)))

## ----eval=FALSE----------------------------------------------------------
#  ##library(diveMove)
#  fname <- system.file(file.path("data", "sealLocs.csv"),
#                         package="diveMove")
#  if (file.exists(fname)) {
#      dat <- read.table(fname, sep=";", header = TRUE, stringsAsFactors = FALSE)
#      dat$class <- ordered(dat$class, c("Z", "B", "A", "0", "1", "2", "3"))
#      dat$time <- as.POSIXct(strptime(dat$time, "%Y-%m-%d %H:%M:%S"), tz = "GMT")
#  }
#  
#  dat1 <- subset(dat, !is.na(dat$lon))
#  coordinates(dat1) <- c("lon", "lat")
#  sl <- trip(dat1, c("time", "id"))
#  
#  ## merge with data from argosfilter
#  #library(argosfilter)
#  #data(seal)
#  ##seal$id <- "ringy2"
#  ##seal[["time"]] <- seal$dtime
#  ##seal$dtime <- NULL
#  ## reconstruct the Argos labels
#  ##seal[["class"]] <- ordered(levels(dat$class)[factor(seal$lc, sort(unique(seal$lc)))], levels(dat$class))
#  ##seal$lc <- NULL
#  
#  
#  ## also adehabitatLT, crawl
#  
#  ## what are we supposed to do with duplicated times?
#  ##  which(!diff(seal$time) > 0)
#  ##[1]   17  116  122 1008 1158 1231 1293 1300
#  ##plot(seal[which(!diff(seal$time) > 0),c("lon", "lat") ])
#  ##points(seal[1 + which(!diff(seal$time) > 0),c("lon", "lat") ], col = "red")
#  
#  ##seal <- seal[!duplicated(seal$time), ]
#  
#  ## drop missing data and combine
#  ##dat <- rbind(dat[!is.na(dat$lon), ], seal[,names(dat)])
#  ##coordinates(dat) <- c("lon", "lat")
#  ##proj4string(dat) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#  
#  ##tr <- trip(dat, c("time", "id"))
#  

