##load data
data(pred.cv.mesa)

##compute long term averages of predictions and observations
pred.lta <- computeLTA(pred.cv.mesa)

##we can now compare observed and predicted averages at each site
plot(pred.lta[,"obs"], pred.lta[,"EX.mu"], pch=1,
     xlim=range(pred.lta), ylim=range(pred.lta),
     xlab="obs", ylab="predictions")
##for the different model components
points(pred.lta[,"obs"], pred.lta[,"EX.mu.beta"], pch=3, col=2)
points(pred.lta[,"obs"], pred.lta[,"EX"], pch=4, col=3)
abline(0,1)

##we could also try computaitons on the original scale
pred.lta <- computeLTA(pred.cv.mesa, exp)

##compare observed and predicted averages
plot(pred.lta[,"obs"], pred.lta[,"EX.mu"], pch=1,
     xlim=range(pred.lta), ylim=range(pred.lta),
     xlab="obs", ylab="predictions")
points(pred.lta[,"obs"], pred.lta[,"EX.mu.beta"], pch=3, col=2)
points(pred.lta[,"obs"], pred.lta[,"EX"], pch=4, col=3)
abline(0,1)
