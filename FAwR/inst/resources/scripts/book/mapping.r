### R code from vignette source 'mapping.rnw'

###################################################
### code chunk number 1: Setup
###################################################
rm(list=ls())

options(repos="http://cran.r-project.org")
options(width=65)

if(!require(class, quietly=TRUE)) install.packages("class")
if(!require(cramer, quietly=TRUE)) install.packages("cramer")
if(!require(fields, quietly=TRUE)) install.packages("fields")
if(!require(gstat, quietly=TRUE)) install.packages("gstat")
if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(lattice, quietly=TRUE)) install.packages("lattice")
if(!require(gpclib, quietly=TRUE)) install.packages("gpclib")
if(!require(maptools, quietly=TRUE)) install.packages("maptools")
if(!require(MASS, quietly=TRUE)) install.packages("MASS")
if(!require(norm, quietly=TRUE)) install.packages("norm")
if(!require(spatial, quietly=TRUE)) install.packages("spatial")
if(!require(spatstat, quietly=TRUE)) install.packages("spatstat")
if(!require(ggplot2, quietly=TRUE)) install.packages("ggplot2")

lattice.options(default.theme = canonical.theme(color = FALSE))

## super hack to keep moving
## sorry andrew, it's the best i could come up with for today
## It'll do just fine, Jeff.
if(getRversion() != "2.9.2") {
  gpclibPermit()
}

##source("../scripts/functions.R")
##options(width=65)                    
#memory.limit(size = 1024*3)



###################################################
### code chunk number 2: mapping.rnw:153-216
###################################################

source("../../scripts/functions.R")
stands <- readShapePoly("../../data/stands.shp", 
                        verbose = FALSE)

## stands.frame <- as.data.frame(stands)
## stands <- stands.frame

## Jeff --- these need to be seen.

## This performs a simple edit, which I'll move out of the book. 
## It should be part of the base data (NA's) rather than zeros (0)
## and I should go ahead and make an external process (without needing the book)
## to generate a unified dataset. 

## Thanks!

## What exactly does the first one do?

## blks <- as.numeric(1:length(stands$STANDID))[is.na(stands$STANDID)]
## stands$STANDID[blks] <- NA

## Could the subsequent ones be like

# stands$TAGE[stands$TAGE == 0] <- NA

# ?  And maybe in a loop?

## blks <- as.numeric(1:length(stands$TAGE))[stands$TAGE == 0]
## stands$TAGE[blks] <- NA

## blks <- as.numeric(1:length(stands$BHAGE))[stands$BHAGE == 0]
## stands$BHAGE[blks] <- NA

## blks <- as.numeric(1:length(stands$DF_SITE))[stands$DF_SITE == 0]
## stands$DF_SITE[blks] <- NA

## blks <- as.numeric(1:length(stands$TOTHT))[stands$TOTHT == 0]
## stands$TOTHT[blks] <- NA

## blks <- as.numeric(1:length(stands$CUBVOL_AC))[stands$CUBVOL_AC == 0]
## stands$CUBVOL_AC[blks] <- NA

## blks <- as.numeric(1:length(stands$TPA))[stands$TPA == 0]
## stands$TPA[blks] <- NA

## blks <- as.numeric(1:length(stands$QMD))[stands$QMD == 0]
## stands$QMD[blks] <- NA

## blks <- as.numeric(1:length(stands$BA))[stands$BA == 0]
## stands$BA[blks] <- NA


## write it out to the disk and make sure you can make ch2 also

## write.polylistShape(polylist, df, file, factor2char = TRUE,
##  strictFilename=FALSE, force = TRUE, max_nchar=254)

##writePolyShape(stands, fn, factor2char = TRUE, max_nchar=254)

##writePolyShape(stands, "stands2", factor2char = TRUE, max_nchar=254)




###################################################
### code chunk number 3: mapping.rnw:219-222
###################################################

show.cols.with.na(as.data.frame(stands))



###################################################
### code chunk number 4: mapping.rnw:237-241
###################################################

stands.frame <- as.data.frame(stands)
stands.frame$MISSING <- rowSums(is.na(stands.frame))



###################################################
### code chunk number 5: mapping.rnw:249-253
###################################################

brks <- sort(unique(stands.frame$MISSING))
colors <- gray(length(brks):1 / (length(brks))) 



###################################################
### code chunk number 6: fig-mcd-attributes
###################################################

plot(stands, 
     col = colors[findInterval(stands.frame$MISSING, 
                               brks, all.inside = TRUE)], 
     forcefill = FALSE,
     axes = TRUE)
title(main="Attribute Missingness by Polygon")
legend(1280000, 365000, brks, fill = colors, cex = 0.7, 
       ncol = 3, title = "# of Missing Attributes")



###################################################
### code chunk number 7: fig-mcd-attributes
###################################################

plot(stands, 
     col = colors[findInterval(stands.frame$MISSING, 
                               brks, all.inside = TRUE)], 
     forcefill = FALSE,
     axes = TRUE)
title(main="Attribute Missingness by Polygon")
legend(1280000, 365000, brks, fill = colors, cex = 0.7, 
       ncol = 3, title = "# of Missing Attributes")



###################################################
### code chunk number 8: mapping.rnw:358-361
###################################################

nap <- na.pattern(stands.frame)



###################################################
### code chunk number 9: mapping.rnw:369-372
###################################################

as.matrix(nap)



###################################################
### code chunk number 10: mapping.rnw:397-460
###################################################

naplot2 <- function(obj, which=c('all','na per var','na per obs','mean na',
                                'na per var vs mean na'),
                   ...)
{
  which <- match.arg(which)
  tab <- table(obj$na.per.obs)
  na.per.var <- diag(obj$sim)
  names(na.per.var) <- dimnames(obj$sim)[[2]]
  mean.na <- obj$mean.na

  if(which %in% c('all','na per var'))
    dotchart(sort(na.per.var), xlab='Fraction of NAs', 
             main='Fraction of NAs in each Variable', ...)

  if(which %in% c('all','na per obs'))
    dotchart2(tab, auxdata=tab, reset.par=TRUE,
              xlab='Frequency', 
              main='Number of Missing Variables Per Observation', ...)

  if(which %in% c('all','mean na'))
    dotchart(sort(mean.na), 
             xlab='Mean Number of NAs',
             main='Mean Number of Other Variables Missing for\nObservations where Indicated Variable is NA',
             ...)

  if(which %in% c('all','na per var vs mean na'))
    {
      if(.R.)
        {
          xpd <- par('xpd')
          par(xpd=NA)
          on.exit(par(xpd=xpd))
        }

      plot(na.per.var, mean.na, xlab='Fraction of NAs for Single Variable',
           ylab='Mean # Other Variables Missing', type='p', ... )
      usr <- par('usr')
      eps <- .015*diff(usr[1:2]);
      epsy <- .015*diff(usr[3:4])
    
      s <- (1:length(na.per.var))[!is.na(mean.na)]
      taken.care.of <- NULL
      for(i in s)
        {
          if(i %in% taken.care.of)
            next

          w <- s[s > i & abs(na.per.var[s]-na.per.var[i]) < eps &
                 abs(mean.na[s]-mean.na[i]) < epsy]
          if(any(w))
            {
              taken.care.of <- c(taken.care.of, w)
              text(na.per.var[i]+eps, mean.na[i],
                   paste(names(na.per.var[c(i,w)]),collapse='\n'),adj=0)
            }
          else text(na.per.var[i]+eps, mean.na[i], names(na.per.var)[i], adj=0)
        }
    }
  
  invisible(tab)
}



###################################################
### code chunk number 11: mapping.rnw:464-467 (eval = FALSE)
###################################################
## par(mfcol=c(2,2), cex=0.7)
## nac <- naclus(stands.frame)
## naplot(nac)


###################################################
### code chunk number 12: fig-nap-plots
###################################################

par(mfcol=c(2,2), cex=0.7)

nac <- naclus(stands.frame)
naplot(nac, which=c('na per var'))
naplot(nac, which=c('na per obs'))
naplot(nac, which=c('mean na'))
naplot2(nac, which=c('na per var vs mean na'), xlim=c(0,0.45))



###################################################
### code chunk number 13: fig-nap-plots
###################################################

par(mfcol=c(2,2), cex=0.7)

nac <- naclus(stands.frame)
naplot(nac, which=c('na per var'))
naplot(nac, which=c('na per obs'))
naplot(nac, which=c('mean na'))
naplot2(nac, which=c('na per var vs mean na'), xlim=c(0,0.45))



###################################################
### code chunk number 14: mapping.rnw:532-537
###################################################

hq <- c(1288538.5625,373896.78125)
centers <- coordinates(stands)
stand.dist <- rdist(t(c(1288538.5625,373896.78125)), centers) 



###################################################
### code chunk number 15: mapping.rnw:543-546
###################################################

site.na <- as.numeric(is.na(stands.frame$DF_SITE))



###################################################
### code chunk number 16: mapping.rnw:551-556
###################################################

gfit.site <- glm(site.na ~ c(stand.dist), 
                 family = "binomial")
summary(gfit.site)



###################################################
### code chunk number 17: mapping.rnw:734-744
###################################################

centers <- coordinates(stands)

stands.frame <- as.data.frame(stands)
stands.frame$x.ctr <- centers[,1]
stands.frame$y.ctr <- centers[,2]
stands.frame$MISSING <- rowSums(is.na(stands.frame))
training.stands <- subset(stands.frame, MISSING == 0)
nrow(training.stands)



###################################################
### code chunk number 18: mapping.rnw:750-754
###################################################

target.stands <- subset(stands.frame, MISSING > 0)
nrow(target.stands)



###################################################
### code chunk number 19: mapping.rnw:771-776
###################################################

cls <- c("x.ctr","y.ctr")
training <- training.stands[,cls]
target <- target.stands[,cls]



###################################################
### code chunk number 20: mapping.rnw:787-791
###################################################

head(rownames(training))
head(rownames(target))



###################################################
### code chunk number 21: mapping.rnw:808-813
###################################################

cl <- factor(rownames(training))

knn.res <- knn(training, target, cl, k = 1)



###################################################
### code chunk number 22: mapping.rnw:834-837
###################################################

missing.rows <- as.character(rownames(target))



###################################################
### code chunk number 23: mapping.rnw:845-848
###################################################

replacement.rows <- as.character(knn.res)



###################################################
### code chunk number 24: mapping.rnw:860-863
###################################################

stands.knn <- stands.frame



###################################################
### code chunk number 25: mapping.rnw:880-886
###################################################

stands.knn[missing.rows,c(4,5:13)] <- 
  stands.frame[replacement.rows,c(4,5:13)]

show.cols.with.na(stands.knn)



###################################################
### code chunk number 26: mapping.rnw:894-899
###################################################

stands.knn$replaced.by <- NA

stands.knn[missing.rows,]$replaced.by <- replacement.rows



###################################################
### code chunk number 27: mapping.rnw:905-910
###################################################

head(as.data.frame(stands.knn))

show.cols.with.na(stands.knn)



###################################################
### code chunk number 28: mapping.rnw:920-923
###################################################

stands.knn[c(1,3),]



###################################################
### code chunk number 29: mapping.rnw:935-938
###################################################

show.cols.with.na(stands.knn)



###################################################
### code chunk number 30: mapping.rnw:1071-1076
###################################################

cls <- c(5:12)
sd <- as.matrix(as.data.frame(stands[,cls]))
psd <- prelim.norm(sd)   



###################################################
### code chunk number 31: mapping.rnw:1090-1093
###################################################

thetahat <- em.norm(psd, showits=FALSE)  #compute mle



###################################################
### code chunk number 32: mapping.rnw:1113-1116
###################################################

em.params <- getparam.norm(psd, thetahat, corr=TRUE)



###################################################
### code chunk number 33: mapping.rnw:1122-1126
###################################################

names(em.params)
em.params$mu



###################################################
### code chunk number 34: mapping.rnw:1131-1134
###################################################

rngseed(1234567)   #set random number generator seed



###################################################
### code chunk number 35: mapping.rnw:1140-1143
###################################################

stands.em <- as.data.frame(stands)



###################################################
### code chunk number 36: mapping.rnw:1149-1152
###################################################

stands.em[,cls] <- imp.norm(psd, thetahat, sd) 



###################################################
### code chunk number 37: mapping.rnw:1157-1160
###################################################

show.cols.with.na(stands.em)



###################################################
### code chunk number 38: mapping.rnw:1233-1246
###################################################

knn.frame <- data.frame(site = stands.knn$DF_SITE, 
                        method = "KNN")

em.frame <- data.frame(site = stands.em$DF_SITE, 
                       method = "EM")

obs.frame <- 
  data.frame(site = subset(stands.frame, 
                           MISSING == 0)$DF_SITE, 
             method = "OBS")




###################################################
### code chunk number 39: mapping.rnw:1252-1255
###################################################

site.frame <- rbind(knn.frame, em.frame, obs.frame)



###################################################
### code chunk number 40: fig-imputed-data-histos
###################################################

histogram(~ as.numeric(site ) | method, data = site.frame,
          xlab = "Site Index (feet)", 
          type = "density", 
          main = "Stand Site Index Frequency Distributions",
          breaks=30,
          layout=c(3,1),
          index.cond=list(c(3,1,2)),
          panel = function(x, ...) {
            panel.histogram(x, ...)
            panel.mathdensity(dmath = dnorm, 
                              col = "black",
                              args = list(mean=mean(x),
                                sd=sd(x)))
          } 
         ) 



###################################################
### code chunk number 41: fig-imputed-data-histos
###################################################
print(

histogram(~ as.numeric(site ) | method, data = site.frame,
          xlab = "Site Index (feet)", 
          type = "density", 
          main = "Stand Site Index Frequency Distributions",
          breaks=30,
          layout=c(3,1),
          index.cond=list(c(3,1,2)),
          panel = function(x, ...) {
            panel.histogram(x, ...)
            panel.mathdensity(dmath = dnorm, 
                              col = "black",
                              args = list(mean=mean(x),
                                sd=sd(x)))
          } 
         ) 

     )


###################################################
### code chunk number 42: mapping.rnw:1515-1518
###################################################

final.plots <- read.csv("../../data/final-plots.csv" )



###################################################
### code chunk number 43: fig-site-plots-hist
###################################################

opar <- par(mfcol = c(1,2), las = 1, cex.axis = 0.70)
site.plots <- subset(final.plots, !is.na(site))
hist(site.plots$site, xlim = c(0,220), ylim = c(0,0.025), 
     freq = FALSE, main = "Histogram", xlab = "Site Index", 
     breaks = 30)
qqnorm(site.plots$site, pch=46)
qqline(site.plots$site, lty=1, col="grey", pch=46)


###################################################
### code chunk number 44: fig-site-plots-hist
###################################################

opar <- par(mfcol = c(1,2), las = 1, cex.axis = 0.70)
site.plots <- subset(final.plots, !is.na(site))
hist(site.plots$site, xlim = c(0,220), ylim = c(0,0.025), 
     freq = FALSE, main = "Histogram", xlab = "Site Index", 
     breaks = 30)
qqnorm(site.plots$site, pch=46)
qqline(site.plots$site, lty=1, col="grey", pch=46)


###################################################
### code chunk number 45: mapping.rnw:1843-1846
###################################################

site.plots <- subset(final.plots, !is.na(site))



###################################################
### code chunk number 46: mapping.rnw:1852-1860
###################################################
site.var <- gstat::variogram(site ~ 1,
                             locations = ~x+y,
                             data = site.plots,
                             width = 50,
                             cutoff = 3000)
site.model <- 
  fit.variogram(site.var, 
                vgm(1000, "Exp", 1000, nugget = 0))


###################################################
### code chunk number 47: mapping.rnw:1870-1873
###################################################

site.model



###################################################
### code chunk number 48: mapping.rnw:1879-1882
###################################################

attributes(site.model)$SSErr



###################################################
### code chunk number 49: fig-site-variogram
###################################################

plot(site.var, 
     model = site.model, 
     main="Site Index Semi-variogram", 
     ylab = "Semi-variance", xlab = "Distance (m)",
     ylim = c(100,310), xlim = c(0,3000))


###################################################
### code chunk number 50: fig-site-variogram
###################################################
print(

plot(site.var, 
     model = site.model, 
     main="Site Index Semi-variogram", 
     ylab = "Semi-variance", xlab = "Distance (m)",
     ylim = c(100,310), xlim = c(0,3000))
     )


###################################################
### code chunk number 51: mapping.rnw:1915-1918
###################################################

sum(site.model[,2])



###################################################
### code chunk number 52: fig-width-cutoff-variograms
###################################################

opar <- par(mfrow=c(2,3), las=1, cex.axis=0.8)

## do the different widths
site.var.10.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=10, cutoff=5000)
site.var.20.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=50, cutoff=5000)
site.var.50.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=100, cutoff=5000)

plot(site.var.10.5000$gamma ~ site.var.10.5000$dist, main="Width=10, Cutoff=5000", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.20.5000$gamma ~ site.var.20.5000$dist, main="Width=50, Cutoff=5000", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.50.5000$gamma ~ site.var.50.5000$dist, main="Width=100, Cutoff=5000", ylab="Semi-variance", xlab="Distance (m)")

## do the different cutoffs
site.var.10.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=10, cutoff=500)
site.var.20.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=50, cutoff=1000)
site.var.50.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=100, cutoff=10000)

plot(site.var.10.5000$gamma ~ site.var.10.5000$dist, main="Width=10, Cutoff=500", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.20.5000$gamma ~ site.var.20.5000$dist, main="Width=50, Cutoff=1000", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.50.5000$gamma ~ site.var.50.5000$dist, main="Width=100, Cutoff=10000", ylab="Semi-variance", xlab="Distance (m)")

par(opar)



###################################################
### code chunk number 53: fig-width-cutoff-variograms
###################################################

opar <- par(mfrow=c(2,3), las=1, cex.axis=0.8)

## do the different widths
site.var.10.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=10, cutoff=5000)
site.var.20.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=50, cutoff=5000)
site.var.50.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=100, cutoff=5000)

plot(site.var.10.5000$gamma ~ site.var.10.5000$dist, main="Width=10, Cutoff=5000", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.20.5000$gamma ~ site.var.20.5000$dist, main="Width=50, Cutoff=5000", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.50.5000$gamma ~ site.var.50.5000$dist, main="Width=100, Cutoff=5000", ylab="Semi-variance", xlab="Distance (m)")

## do the different cutoffs
site.var.10.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=10, cutoff=500)
site.var.20.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=50, cutoff=1000)
site.var.50.5000 <- gstat::variogram(site~1,locations=~x+y,data=site.plots,width=100, cutoff=10000)

plot(site.var.10.5000$gamma ~ site.var.10.5000$dist, main="Width=10, Cutoff=500", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.20.5000$gamma ~ site.var.20.5000$dist, main="Width=50, Cutoff=1000", ylab="Semi-variance", xlab="Distance (m)")
plot(site.var.50.5000$gamma ~ site.var.50.5000$dist, main="Width=100, Cutoff=10000", ylab="Semi-variance", xlab="Distance (m)")

par(opar)



###################################################
### code chunk number 54: mapping.rnw:1987-1990
###################################################

bdry2 <- readShapePoly("../../data/boundary.shp")



###################################################
### code chunk number 55: mapping.rnw:1997-2001
###################################################

dum <- fortify(bdry2, region = "FORBNDRY_")
summary(dum)



###################################################
### code chunk number 56: mapping.rnw:2019-2025
###################################################

dum2 <- dum[dum$id == 3,]
outer <- dum2[dum2$piece == 1,]
inner <- dum2[dum2$piece == 2,]




###################################################
### code chunk number 57: mapping.rnw:2034-2041
###################################################
bdry.poly <- vector(2, mode="list")
bdry.poly[[1]] <- 
  list(x = outer[(nrow(outer)-1):1,]$long, 
       y = outer[(nrow(outer)-1):1,]$lat)
bdry.poly[[2]] <- 
  list(x = inner[(nrow(inner)-1):1,]$long, 
       y = inner[(nrow(inner)-1):1,]$lat)


###################################################
### code chunk number 58: mapping.rnw:2047-2051
###################################################

bdry.owin <- owin(poly=bdry.poly)
grid <- gridcentres(bdry.owin, 200, 200) 



###################################################
### code chunk number 59: mapping.rnw:2065-2069
###################################################

ok <- inside.owin(grid$x, grid$y, bdry.owin)
pred.grid <- data.frame(x=grid$x[ok], y=grid$y[ok])



###################################################
### code chunk number 60: mapping.rnw:2075-2083
###################################################

sample.site.plots <- sample(1:nrow(site.plots), 1000 )
site.pred <- krige(formula = site ~ 1,
                  locations = ~ x + y,
                  data = site.plots[sample.site.plots,],
                  model = site.model,
                  newdata = pred.grid)



###################################################
### code chunk number 61: mapping.rnw:2091-2094
###################################################

head(site.pred)



###################################################
### code chunk number 62: fig-site-index-histos
###################################################

opar <- par(mfcol = c(1,2), las = 1, cex.axis = 0.80)
hist(site.pred$var1.pred, breaks=30, 
     main="Predicted Site Index", 
     freq=FALSE, ylim=c(0,0.08), xlim=c(80,150))
qqnorm(site.pred$var1.pred, pch=46)
qqline(site.pred$var1.pred, lty=1, col="darkgrey", pch=46)


###################################################
### code chunk number 63: fig-site-index-histos
###################################################

opar <- par(mfcol = c(1,2), las = 1, cex.axis = 0.80)
hist(site.pred$var1.pred, breaks=30, 
     main="Predicted Site Index", 
     freq=FALSE, ylim=c(0,0.08), xlim=c(80,150))
qqnorm(site.pred$var1.pred, pch=46)
qqline(site.pred$var1.pred, lty=1, col="darkgrey", pch=46)


###################################################
### code chunk number 64: fig-site-index-pred
###################################################
levelplot(var1.pred ~ x + y,
          data = site.pred,
          col.regions = terrain.colors(80),
          main = "Predicted Site Index",
          xlab = "Easting",
          ylab = "Northing",
          panel = function(...) {
            panel.levelplot(...)
            lpoints(site.plots[sample.site.plots,]$x, 
                    site.plots[sample.site.plots,]$y, 
                    col="black", cex=0.25, pch=19)
          }
         )


###################################################
### code chunk number 65: fig-site-index-pred
###################################################
print(
levelplot(var1.pred ~ x + y,
          data = site.pred,
          col.regions = terrain.colors(80),
          main = "Predicted Site Index",
          xlab = "Easting",
          ylab = "Northing",
          panel = function(...) {
            panel.levelplot(...)
            lpoints(site.plots[sample.site.plots,]$x, 
                    site.plots[sample.site.plots,]$y, 
                    col="black", cex=0.25, pch=19)
          }
         )
     )


###################################################
### code chunk number 66: fig-site-index-variance
###################################################
print(
levelplot(var1.var ~ x+y,
          data = site.pred,
          col.regions = terrain.colors(80),
          main = "Estimated Variance of Site Index Prediction",
          xlab = "Easting",
          ylab = "Northing",
          panel = function(...) {
            panel.levelplot(...)
            lpoints(site.plots[sample.site.plots,]$x, 
                    site.plots[sample.site.plots,]$y, 
                    col="black", cex=0.25, pch=19)
          }
         )
     )


###################################################
### code chunk number 67: mapping.rnw:2521-2523
###################################################
system("rm -fr package-Ch4")
package.skeleton(name = "package-Ch4")


