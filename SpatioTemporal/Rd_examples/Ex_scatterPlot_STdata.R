################################
## Example for STdata/STmodel ##
################################
##load data
data(mesa.model)

par(mfrow=c(2,2))
##plot observations as a function of longitude for an STmodel object
scatterPlot(mesa.model, covar="long")

##as a function of the first temporal trend, subset to only AQS sites
##and fit for each location
scatterPlot(mesa.model, trend=1, col=c(1:25,1), pch=19, cex=.1,
            group=mesa.model$obs$ID, lty=c(rep(2,25),1),
            subset=with(mesa.model$locations, ID[type=="AQS"]))

##if plotting against the distance to coast, we might have to change the
##smooting.
suppressWarnings( scatterPlot(mesa.model, covar="km.to.coast") )
##better
scatterPlot(mesa.model, covar="km.to.coast", col=c(NA,2), add=TRUE,
            smooth.args=list(span=4/5,degree=2))

##Lets group data by season
##First create a vector dividing data into four seasons
I.season <- as.factor(as.POSIXlt(mesa.model$obs$date)$mon+1)
levels(I.season) <- c(rep("Winter",2), rep("Spring",3), 
                      rep("Summer",3), rep("Fall",3), "Winter") 
scatterPlot(mesa.model, covar="log10.m.to.a1", col=c(2:5,1),
            group=I.season)
legend("bottomleft", c(levels(I.season),"All"), col=c(2:5,1), pch=1)


###############################
## Example for predCVSTmodel ##
###############################
##load data
data(pred.cv.mesa)

##simple case of residuals against temporal trends
par(mfrow=c(2,1))
scatterPlot(pred.cv.mesa, trend=1, STdata=mesa.model, type="res")

##colour coded by season
I.season <- as.factor(as.POSIXlt(pred.cv.mesa$pred.obs$date)$mon+1)
levels(I.season) <- c(rep("Winter",2), rep("Spring",3), 
                      rep("Summer",3), rep("Fall",3), "Winter") 

scatterPlot(pred.cv.mesa, trend=1, STdata=mesa.model, type="res",
            group=I.season, col=c(2:5,1), lty=c(1,1,1,1,2),
            smooth.args=list(span=.1,degree=2))
            
##or as function of covariates
par(mfcol=c(2,2))
scatterPlot(pred.cv.mesa, , type="res", covar="log10.m.to.a1",
            STdata=mesa.model, group=I.season, col=c(2:5,1))
scatterPlot(pred.cv.mesa, type="res", covar="km.to.coast",
            STdata=mesa.model, group=I.season, col=c(2:5,1),
            smooth.args=list(span=4/5,degree=1))

##let's compare to the original observations
scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", STdata=mesa.model,
            group=I.season, col=c(2:5,1), type="obs")
scatterPlot(pred.cv.mesa, covar="km.to.coast", STdata=mesa.model,
            group=I.season, col=c(2:5,1), type="obs",
            smooth.args=list(span=4/5,degree=1))
