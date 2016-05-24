## ----message=FALSE-------------------------------------------------------
library(crawl)

## ------------------------------------------------------------------------
data("northernFurSeal")
head(northernFurSeal)

## ------------------------------------------------------------------------
northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class,
                                          levels=c("3", "2", "1","0","A"))

## ---- message=FALSE------------------------------------------------------
library(sp)
library(rgdal)
coordinates(northernFurSeal) = ~longitude+latitude

## ------------------------------------------------------------------------
proj4string(northernFurSeal) <- CRS("+proj=longlat")

## ----message=FALSE-------------------------------------------------------
northernFurSeal <- spTransform(northernFurSeal, 
                               CRS(paste("+proj=aea +lat_1=30 +lat_2=70",
                                         "+lat_0=52 +lon_0=-170 +x_0=0 +y_0=0",
                                         "+ellps=GRS80 +datum=NAD83",
                                         "+units=m +no_defs"))
)

## ----message=FALSE-------------------------------------------------------
initial = list(a=c(coordinates(northernFurSeal)[1,1],0,
                   coordinates(northernFurSeal)[1,2],0),
               P=diag(c(10000^2,54000^2,10000^2,5400^2)))

## ----message=FALSE-------------------------------------------------------
fixPar = c(log(250), log(500), log(1500), rep(NA,3), NA)

## ----message=FALSE-------------------------------------------------------
displayPar(mov.model=~1,
           err.model=list(x=~Argos_loc_class-1),
           data=northernFurSeal,
           fixPar=fixPar)

## ---- message=FALSE------------------------------------------------------
constr=list(lower=c(rep(log(1500),2), rep(-Inf,2)),
            upper=rep(Inf,4))

## ---- message=FALSE------------------------------------------------------
ln.prior = function(theta){-abs(theta[4]-3)/0.5}

## ----message=FALSE-------------------------------------------------------
set.seed(123)
fit1 <- crwMLE(mov.model=~1, 
               err.model=list(x=~Argos_loc_class-1),
               data=northernFurSeal, 
               Time.name="Time",
               initial.state=initial,
               fixPar=fixPar, 
               constr=constr, 
               prior=ln.prior,
               control=list(maxit=30, trace=0,REPORT=1),
               initialSANN=list(maxit=200, trace=0, REPORT=1))

## ------------------------------------------------------------------------
fit1

## ----message=FALSE-------------------------------------------------------
predTime <- seq(ceiling(min(northernFurSeal$Time)), 
                floor(max(northernFurSeal$Time)), 1)

## ----message=FALSE-------------------------------------------------------
predObj <- crwPredict(object.crwFit=fit1, 
                      predTime, 
                      speedEst=TRUE, 
                      flat=TRUE)

## ----message=FALSE-------------------------------------------------------
crwPredictPlot(predObj, "map")

## ----message=FALSE-------------------------------------------------------
set.seed(123)
simObj <- crwSimulator(fit1, 
                       predTime, 
                       method="IS", 
                       parIS=100, 
                       df=5, 
                       scale=18/20)

## ----message=FALSE-------------------------------------------------------
w <- simObj$thetaSampList[[1]][,1]
hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')

## ----message=FALSE-------------------------------------------------------
round(100/(1+(sd(w)/mean(w))^2))

## ------------------------------------------------------------------------
my.colors <-colorRampPalette(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))

## ------------------------------------------------------------------------
iter <- 20
cols <- my.colors(iter)

## ----message=FALSE-------------------------------------------------------
crwPredictPlot(predObj, 'map')
for(i in 1:iter){
  samp <- crwPostIS(simObj)
  lines(samp$alpha.sim[,'mu.x'], samp$alpha.sim[,'mu.y'],col=cols[i]) 
}

## ------------------------------------------------------------------------
library(crawl)
library(sp)
library(rgdal)
library(ggplot2)
data("harborSeal")
head(harborSeal)

## ------------------------------------------------------------------------
harborSeal$Argos_loc_class = factor(harborSeal$Argos_loc_class, levels=c("3","2","1","0","A","B"))

## ------------------------------------------------------------------------
toProj = harborSeal[!is.na(harborSeal$latitude),c("Time","latitude","longitude")]

## ---- message=FALSE------------------------------------------------------
coordinates(toProj) = ~longitude+latitude

## ------------------------------------------------------------------------
proj4string(toProj) <- CRS("+proj=longlat")
toProj <- spTransform(toProj, CRS("+init=epsg:3338"))

## ------------------------------------------------------------------------
toProj = as.data.frame(toProj)
colnames(toProj)[2:3] = c("x","y")
harborSeal = merge(toProj, harborSeal, by="Time", all=TRUE)
harborSeal = harborSeal[order(harborSeal$Time),]

## ------------------------------------------------------------------------
initial = list(a=c(harborSeal$x[1],0,harborSeal$y[1],0),
               P=diag(c(10000^2,5400^2,10000^2,5400^2)))

## ------------------------------------------------------------------------
fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)
displayPar( mov.model=~1, err.model=list(x=~Argos_loc_class-1),data=harborSeal,activity=~I(1-DryTime),fixPar=fixPar)

## ------------------------------------------------------------------------
constr=list(lower=c(rep(log(1500),3), rep(-Inf,2)),
            upper=rep(Inf,5))

## ----message=FALSE-------------------------------------------------------
set.seed(123)
fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal, coord=c("x","y"), Time.name="Time", 
  initial.state=initial, fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1)
)

## ------------------------------------------------------------------------
print(fit1)

## ------------------------------------------------------------------------
pred1 = crwPredict(fit1, predTime=NULL, flat=TRUE)

## ------------------------------------------------------------------------
p1=ggplot(aes(x=mu.x, y=mu.y), data=pred1) + geom_path(col="red") + geom_point(aes(x=x, y=y), col="blue") + coord_fixed()

p2=ggplot(aes(x=Time, y=mu.x), data=pred1) + geom_ribbon(aes(ymin=mu.x-2*se.mu.x,ymax=mu.x+2*se.mu.x),fill="green", alpha=0.5)  + geom_path(col="red") + geom_point(aes(x=Time, y=x), col="blue", size=1)

p3=ggplot(aes(x=Time, y=mu.y), data=pred1) + 
  geom_ribbon(aes(ymin=mu.y-2*se.mu.y,ymax=mu.y+2*se.mu.y), fill="green", alpha=0.5)  + 
  geom_path(col="red") + geom_point(aes(x=Time, y=y), col="blue", size=1)

suppressWarnings(print(p1))
suppressWarnings(print(p2))
suppressWarnings(print(p3))

## ----message=FALSE-------------------------------------------------------
library(sp)
library(rgdal)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)

## ----message=FALSE-------------------------------------------------------
data("beardedSeals")
beardedSeals

## ----message=FALSE-------------------------------------------------------
beardedSeals %>% 
  group_by(deployid,date_time) %>% 
  filter(n()>1)

## ----message=FALSE-------------------------------------------------------
library(xts)
date_unique <-beardedSeals %>% 
  group_by(deployid) %>%
  do(unique_date = xts::make.time.unique(.$date_time,eps=1)) %>%
  tidyr::unnest(unique_date) %>%
  mutate(unique_posix = as.POSIXct(.$unique_date,origin='1970-01-01 00:00:00',tz='UTC')) %>%
  dplyr::arrange(deployid,unique_posix) %>% 
  dplyr::select(unique_posix)

beardedSeals <- beardedSeals %>% arrange(deployid,date_time) %>%
  bind_cols(date_unique)

## ----message=FALSE-------------------------------------------------------
beardedSeals %>% 
  group_by(deployid,unique_posix) %>% 
  filter(n()>1)

## ----message=FALSE-------------------------------------------------------
beardedSeals <- beardedSeals %>%  
  dplyr::arrange(deployid,unique_posix)

library(doParallel)
library(argosfilter)

split_data <- split(beardedSeals,beardedSeals$deployid)

registerDoParallel(cores=2)
beardedSeals$filtered <- foreach(i = 1:length(split_data), .combine = c) %dopar% {
  argosfilter::sdafilter(
    lat=split_data[[i]]$latitude, 
    lon=split_data[[i]]$longitude, 
    dtime=split_data[[i]]$unique_posix,
    lc=split_data[[i]]$quality, 
    ang=-1,
    vmax=5)
}
stopImplicitCluster()

beardedSeals <- beardedSeals %>% 
  dplyr::filter(., filtered=="not" & !is.na(error_semimajor_axis)) %>%
  arrange(.,deployid,unique_posix)

## ----message=FALSE-------------------------------------------------------
beardedSeals <- as.data.frame(beardedSeals)
coordinates(beardedSeals) = ~longitude+latitude
proj4string(beardedSeals) = CRS("+proj=longlat +datum=WGS84")

beardedSeals <- spTransform(beardedSeals, CRS("+init=epsg:3571"))

## ----message=FALSE-------------------------------------------------------
ids = unique(beardedSeals@data$deployid)      #define seal IDs

registerDoParallel(cores=2)
model_fits <-
  foreach(i = 1:length(ids)) %dopar% {
    id_data = subset(beardedSeals,deployid == ids[i])
    diag_data = model.matrix(
      ~ error_semimajor_axis + error_semiminor_axis + error_ellipse_orientation,
      id_data@data
    )[,-1]
    
    id_data@data = cbind(id_data@data, 
                         crawl::argosDiag2Cov(
                           diag_data[,1], 
                           diag_data[,2], 
                           diag_data[,3]))
    
    init = list(a = c(sp::coordinates(id_data)[1,1],0,
                      sp::coordinates(id_data)[1,2],0),
                P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                           5000 ^ 2, 10 * 3600 ^ 2)))
    
    fit <- crawl::crwMLE(
      mov.model =  ~ 1,
      err.model = list(
        x =  ~ ln.sd.x - 1, 
        y =  ~ ln.sd.y - 1, 
        rho =  ~ error.corr
      ),
      data = id_data,
      Time.name = "unique_posix",
      initial.state = init,
      fixPar = c(1,1,NA,NA),
      theta = c(log(10), 3),
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1)
    )
    fit
  }
stopImplicitCluster()

names(model_fits) <- ids

print(model_fits)

## ----message=FALSE-------------------------------------------------------
registerDoParallel(cores=2)
predData <- foreach(i = 1:length(model_fits), .combine = rbind) %dopar% {
  
  model_fits[[i]]$data$unique_posix <- lubridate::with_tz(
    model_fits[[i]]$data$unique_posix,"GMT")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits[[i]]$data$unique_posix),"hour"),
    lubridate::floor_date(max(model_fits[[i]]$data$unique_posix),"hour"),
    "1 hour")
  tmp = crawl::crwPredict(model_fits[[i]], predTime=predTimes)
}
stopImplicitCluster()

predData$predTimes <- intToPOSIX(predData$TimeNum)

## ----plot-1--------------------------------------------------------------
theme_map = function(base_size=9, base_family="")
{
  require(grid)
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.title.x=element_text(vjust=0),
          axis.title.y=element_text(angle=90, vjust=1.25),
          axis.text.y=element_text(angle=90),
          axis.ticks=element_line(colour="black", size=0.25),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text=element_text(),
          legend.title=element_text(face="bold", hjust=0),
          panel.border=element_rect(fill=NA, colour="black"),
          panel.grid.major=element_line(colour="grey92", size=0.3, linetype=1),
          panel.grid.minor=element_blank(),
          plot.title=element_text(vjust=1),
          strip.background=element_rect(fill="grey90", colour="black", size=0.3),
          strip.text=element_text()
    )
}

p1 <- ggplot(data=predData,aes(x=mu.x,y=mu.y)) + 
  geom_path(aes(colour=deployid)) + xlab("easting (meters)") +
  ylab("northing (meters)") + theme_map()
p1

