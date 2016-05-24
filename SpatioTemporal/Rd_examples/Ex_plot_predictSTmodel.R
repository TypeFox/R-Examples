#######################################
## plot predictions for a given site ##
#######################################
##load data
data(mesa.model)
##load predictions
data(pred.mesa.model)

par(mfrow=c(2,1))
plot(pred.mesa.model)
##different site with data and prediction variances
plot(pred.mesa.model, STmodel=mesa.model, ID="L001",
     pred.var=TRUE)

##compare the different contributions to the predictions
plot(pred.mesa.model)
plot(pred.mesa.model, pred.type="EX.mu", col="red", add=TRUE)
plot(pred.mesa.model, pred.type="EX.mu.beta", col="green", add=TRUE)

##compare the two confidence and prediction intervalls
plot(pred.mesa.model, ID=3, pred.var=TRUE, col=c(0,0,"darkgrey"))
plot(pred.mesa.model, ID=3, STmodel=mesa.model,
     col=c("black","red","lightgrey"), add=TRUE)

##plot predictions as function of observations
par(mfrow=c(2,2))
plot(pred.mesa.model, y="obs", STmodel=mesa.model, pred.var=TRUE)

##all data, using points and colour coded by site
plot(pred.mesa.model, y="obs", STmodel=mesa.model, ID="all",
     lty=c(NA,1), pch=c(19,NA), col=c("ID", "red", "grey"),
     cex=.25, pred.var=TRUE)

##compare prediction methods, for one site only
plot(pred.mesa.model, y="obs", STmodel=mesa.model,
     lty=c(NA,1), pch=c(19,NA), cex=.25, pred.var=TRUE)
plot(pred.mesa.model, y="obs", STmodel=mesa.model, col="red",
     lty=NA, pch=c(19,NA), cex=.25, pred.type="EX.mu",
     add=TRUE)
plot(pred.mesa.model, y="obs", STmodel=mesa.model, col="green",
     lty=NA, pch=c(19,NA), cex=.25, pred.type="EX.mu.beta",
     add=TRUE)

####################################
## plot CV-pred. for a given site ##
####################################
##load CV-predictions
data(pred.cv.mesa)

par(mfcol=c(3,1),mar=c(2.5,2.5,2,.5))
plot(pred.cv.mesa, ID=1)
plot(pred.cv.mesa, ID=1, pred.type="EX.mu", col="green", add=TRUE)
plot(pred.cv.mesa, ID=1, pred.type="EX.mu.beta", col="blue", add=TRUE)

##different colours
plot(pred.cv.mesa, ID=10, col=c("blue","magenta","light blue"))

##points and lines for the observations
plot(pred.cv.mesa, ID=17, lty=c(1,NA), pch=c(NA,19), cex=.5)

##plot predictions as function of observations
par(mfrow=c(2,2))
plot(pred.cv.mesa, y="obs")

##all data, using points and colour coded by site
plot(pred.cv.mesa, y="obs", ID="all", lty=c(NA,1),
     pch=c(19,NA), cex=.25, col=c("ID", "red", "grey"))

##compare prediction methods, for one site only
plot(pred.cv.mesa, y="obs", lty=c(NA,1), pch=c(19,NA), cex=.25)
plot(pred.cv.mesa, y="obs", col="red", lty=NA, pch=c(19,NA),
     cex=.25, pred.type="EX.mu", add=TRUE)
plot(pred.cv.mesa, y="obs", col="green", lty=NA, pch=c(19,NA),
     cex=.25, pred.type="EX.mu.beta", add=TRUE)
