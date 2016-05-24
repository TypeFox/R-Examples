
################################################################
###################  New York data example #####################
################################################################

# define a pause function
pause <- function() invisible(readline("\n Please Enter to continue: \n"))


# load libraries
library("spTimer");
library("maps");
library("colorspace");

# Read data 
data(NYdata)
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 8)))


# Figure 7
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
map(database="state",regions="new york")
points(coords,pch=19,col=3)
points(coords,pch=1,col=1)
points(pred.coords,pch=3,col=4)
legend(x=-77.5,y=41.5,col=c(3,4),pch=c(19,3),cex=0.8,legend=c("Fitted sites","Validation sites"))


# Fit GP model 
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,data=DataFit, 
        model="GP", coords=~Longitude+Latitude, scale.transform="SQRT",
		spatial.decay=spT.decay(distribution=Gamm(2,1),tuning=0.1))
print(post.gp)
summary(post.gp)


# Spatial prediction for the GP model
set.seed(11)
pred.gp <- predict(post.gp, newdata=DataValPred, newcoords=~Longitude+Latitude)
print(pred.gp)
names(pred.gp)

# model summary
summary(post.gp)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(pred.gp$Median))  


###############################
## For surface plots         ##
## Press Enter:              ##
###############################
pause()

nItr=100
nBurn=50
# Predict on grids
data(NYgrid)
set.seed(11)
post.gp2 <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,data=NYdata, 
        model="GP", coords=~Longitude+Latitude, scale.transform="SQRT",
		spatial.decay=spT.decay(distribution=Gamm(2,1),tuning=0.1))
set.seed(11)
grid.pred <- predict(post.gp2, newdata=NYgrid, newcoords=~Longitude+Latitude)

# predictive plots
library(MBA)
library(fields)
library(maps)

# this function is used to delete values outside NY
fnc.delete.map.XYZ<-function(xyz){
	x<-xyz$x; y<-xyz$y; z<-xyz$z
	xy<-expand.grid(x, y)
	eus<-(map.where(database="state", x=xy[,1], y=xy[,2]))
	dummy<-rep(0, length(xy[,1]))
	eastUS<-NULL
	eastUS<-data.frame(lon=xy[,1],lat=xy[,2],state=eus,dummy=dummy)
	eastUS[!is.na(eastUS[,3]),4]<-1
	eastUS[eastUS[,3]=="pennsylvania" & !is.na(eastUS[,3]),4]<-0
	eastUS[eastUS[,3]=="new jersey" & !is.na(eastUS[,3]),4]<-0
	eastUS[eastUS[,3]=="connecticut" & !is.na(eastUS[,3]),4]<-0
	eastUS[eastUS[,3]=="massachusetts:main" & !is.na(eastUS[,3]),4]<-0
	eastUS[eastUS[,3]=="new hampshire" & !is.na(eastUS[,3]),4]<-0
	eastUS[eastUS[,3]=="vermont" & !is.na(eastUS[,3]),4]<-0
	a <- eastUS[, 4]
	z <- as.vector(xyz$z)
	z[!a] <- NA
	z <- matrix(z, nrow = length(xyz$x))
      xyz$z <- z
      xyz
}
##

coords<-unique(NYdata[,c("Longitude","Latitude")])
grid.coords<-unique(NYgrid[,c("Longitude","Latitude")])
true.val<-matrix(NYdata$o8hrmax,62,28)
grid.val<-matrix(grid.pred$Median,62,dim(grid.coords)[[1]])
grid.sd<-matrix(grid.pred$SD,62,dim(grid.coords)[[1]])

surfplot<-function(day=60, val, ...)
{
    z <- val
	surf<-cbind(grid.coords,z[day,])
	surf<-mba.surf(surf,200,200)$xyz
	surf<-fnc.delete.map.XYZ(xyz=surf)
	#map(database="state",regions="new york")
	image.plot(surf, xlab="Longitude",ylab="Latitude",axes=F, ...)
	contour(surf,nlevels=10,lty=3,add=T)
	map(database="state",regions="new york",add=T)
	axis(1);axis(2)
}

# prediction for day 60
day<-60
surfplot(day, val=grid.val, col = rainbow_hcl(100, start = 200, end = 0))
text(coords,labels=round(true.val[day,],1),cex=0.8,col=1)

# sd for day 60
# Press Enter:
pause()

# sd for day 60
day<-60
surfplot(day, val=grid.sd,col = diverge_hcl(100, h = c(246, 40), c = 96, l = c(65, 90)))
points(coords,pch=19,cex=1,col=2)
points(coords,pch=1,cex=1,col=1)


###############################################################
## To run more code on model and MCMC diagnosis press Enter: ##
###############################################################
pause()

rm(post.gp2);

# More code for summary and plots
summary(post.gp, digits = 4)
summary(post.gp,pack="coda") # mcmc summary statistics using coda package
confint(post.gp)
# MCMC chains
plot(post.gp)
# residual plot
plot(post.gp, residuals=TRUE)
# fitted surface 
library(akima)
plot.spT<-function(x, residuals=FALSE, surface=NULL, time=c(1), a3d=FALSE, 
            points=FALSE, title=TRUE, ...){
   if(is.null(surface) & a3d==FALSE){
      if(as.logical(residuals)==FALSE){
         tmp<-as.mcmc(x)
         plot(tmp, ...)}
      else{
         plot(x$fitted[,1],residuals(x),ylab="Residuals",xlab="Fitted values")
         abline(h=0,lty=2);title("Residuals vs Fitted")
         par(ask=TRUE)
         qqnorm(residuals(x));qqline(residuals(x),lty=2)}} 
   else {
      if(is.null(surface)){
         stop("\n# Error: surface should be defined as 'Mean' or 'SD'. \n")}	
      if(!surface %in% c("Mean","SD")){
         stop("\n# Error: surface only takes 'Mean' or 'SD'. \n")}
      library(akima); library(fields); 
	     z<-array(fitted(x)[,paste(surface)],dim=c(x$T*x$r,x$n))
      xyz<-cbind(x$coords,c(z[time,]))
      xyz<-interp(x=xyz[,1],y=xyz[,2],z=xyz[,3],
           xo=seq(min(xyz[,1]),max(xyz[,1]),length=150),
		        yo=seq(min(xyz[,2]), max(xyz[,2]), length = 150))
      if(a3d==TRUE){
         persp(x=xyz$x,y=xyz$y,z=xyz$z, xlab="x",ylab="y",zlab="z", ...)->res}
	     else{
	        image.plot(xyz, ...) 
	        if(points != FALSE){
	           points(x$coords,pch=16,cex=0.8)}}
	     if(title==TRUE){
	        title(main=paste("Time point: (t=",time,")",sep=""))}}}
contour.spT<-function(x, surface="Mean", time=c(1), ...){
	z<-array(fitted(x)[,paste(surface)],dim=c(x$T*x$r,x$n))
	xyz<-cbind(x$coords,c(z[time,]))
	xyz<-interp(x=xyz[,1], y=xyz[,2], z=xyz[,3], xo=seq(min(xyz[,1]), max(xyz[,1]), length = 150),
		 yo=seq(min(xyz[,2]), max(xyz[,2]), length = 150),linear = TRUE, extrap=FALSE, 
		 duplicate = "error", dupfun = NULL, ncp = NULL)
	contour(xyz, ...) 
}
plot(post.gp, surface="Mean")
# fitted surface 3d plot
plot(post.gp, surface="Mean", a3d=TRUE)
# some other R functions 
coef(post.gp)
formula(post.gp)
terms(post.gp)
head(model.frame(post.gp))
head(model.matrix(post.gp))
# Model selection criteria
post.gp$PMCC 
# MCMC diagnostics using coda
# autocorr diagnostics
autocorr.diag(as.mcmc(post.gp))
# Raftery and Lewis's diagnostic
raftery.diag(post.gp)
# Diagnostics using more than one chain
set.seed(22)
post.gp2 <- NULL
post.gp2 <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,data=DataFit, 
        model="GP", coords=~Longitude+Latitude, scale.transform="SQRT",
        initials=spT.initials(model="GP",phi=1,sig2eta=1),
		spatial.decay=spT.decay(distribution=Gamm(2,1),tuning=0.1))
mcobj<-list(as.mcmc(post.gp),as.mcmc(post.gp2))
mcobj<-as.mcmc.list(mcobj)
# acf plot
acfplot(mcobj)
# Geweke's convergence diagnostic
geweke.diag(mcobj)
# Gelman and Rubin's diagnostic
gelman.diag(mcobj)
gelman.plot(mcobj)


#####################################################
## For Temporal prediction/forecast using GP model ##
## 1. In the unobserved locations                  ## 
## Press Enter:                                    ##
#####################################################
pause()


# Temporal  prediction/forecast for the GP model
# 1. In the unobserved locations
# Read data
data(NYdata)
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day %in% c(30, 31) & Month == 8)))
# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.gp <- predict(post.gp, newdata=DataValFore, newcoords=~Longitude+Latitude, 
        type="temporal", foreStep=2)
print(fore.gp)
names(fore.gp)
# Forecast validations 
spT.validation(DataValFore$o8hrmax,c(fore.gp$Median)) 


#####################################################
## For Temporal prediction/forecast using GP model ##
## 2. In the observed/fitted locations             ##
## Press Enter:								       ##
#####################################################
pause()


# Temporal  prediction/forecast for the GP model
# 2. In the observed/fitted locations
# Read data
s <-c(8,11,12,14,18,21,24,28)
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day %in% c(30, 31) & Month == 8)))


# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.gp <- predict(post.gp, newdata=DataFitFore, newcoords=~Longitude+Latitude, 
           type="temporal", foreStep=2)
print(fore.gp)
names(fore.gp)
# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.gp$Median)) # 


#######################################
## Model data with spacetime classes ##
#######################################
pause()

rm(post.gp); rm(post.gp2); rm(pred.gp); rm(fore.gp); rm(mcobj)

# Create a dataset with spacetime class
library(spTimer)
site<-unique(NYdata[,c("Longitude","Latitude")])
library(spacetime)
row.names(site)<-paste("point",1:nrow(site),sep="")
site <- SpatialPoints(site)
ymd<-as.POSIXct(seq(as.Date("2006-07-01"),as.Date("2006-08-31"),by=1))
# introduce class STFDF
newNYdata<-STFDF(sp=site, time=ymd, data=NYdata) # full lattice
class(newNYdata)

# model with GP
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        data=newNYdata, model="GP", scale.transform="SQRT")
summary(post.gp)


############################################
## Model spatial only data using GP model ##
############################################
pause()


# spatial only data 
# we use meuse data from sp package
library(sp)
data(meuse)
# model with GP
set.seed(11)
post.gp <- spT.Gibbs(formula=zinc ~ sqrt(dist),   
    data=meuse, model="GP", coords=~x+y, nItr=500, nBurn=100,
	spatial.decay=spT.decay(distribution=Gamm(2,1), tuning=0.5),
	distance.method="euclidean",scale.transform="LOG")
summary(post.gp)

plot(post.gp, surface="Mean", title=FALSE)


###################################
## To run more code on AR model: ##
## Press Enter:                  ##
###################################
pause()

rm(post.gp); 

##  AR model
# MCMC via Gibbs
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
          data=DataFit, model="AR", coords=~Longitude+Latitude, 
          scale.transform="SQRT")
print(post.ar)
# Summary and plots
summary(post.ar)
summary(post.ar,pack="coda")
plot(post.ar)
plot(post.ar,residuals=TRUE)
# Model selection criteria
post.ar$PMCC 


# Spatial prediction/interpolation for the AR model
# Define prediction coordinates
set.seed(11)
pred.ar <- predict(post.ar, newdata=DataValPred, newcoords=~Longitude+Latitude)
print(pred.ar)
names(pred.ar)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(pred.ar$Median))  


####################################################
## Temporal  prediction/forecast for the AR model ##
## 1. In the unobserved locations                 ##
## Press Enter:                                   ##
####################################################
pause()


# Temporal  prediction/forecast for the AR model
# 1. In the unobserved locations
# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.ar <- predict(post.ar, newdata=DataValFore, newcoords=~Longitude+Latitude, 
           type="temporal", foreStep=2, predAR=pred.ar)
print(fore.ar)
names(fore.ar)
# Forecast validations 
spT.validation(DataValFore$o8hrmax,c(fore.ar$Median)) 


####################################################
## Temporal  prediction/forecast for the AR model ##
## 2. In the observed/fitted locations            ## 
## Press Enter:                                   ##
####################################################
pause()


# Temporal  prediction/forecast for the AR model
# 2. In the observed/fitted locations
# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.ar <- predict(post.ar, newdata=DataFitFore, newcoords=~Longitude+Latitude, 
           type="temporal", foreStep=2)
print(fore.ar)
names(fore.ar)
# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.ar$Median)) # 


################################# The End ########################################
