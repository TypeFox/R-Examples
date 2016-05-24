##load data
data(mesa.model)

##naive predictions based on either AQS,
pred.aqs <- predictNaive(mesa.model, type="AQS")
##...or only one sites,
pred.1site <- predictNaive(mesa.model, locations="60372005")

##plot the predictions - The two cases that are constant in space
par(mfcol=c(2,1), mar=c(4.5,4.5,1,.5))

##observations as a function of date
plot(mesa.model, "loc.obs", type=as.factor(mesa.model$locations$ID),
     legend.loc=NULL, pch=19, cex=.25)
##Add the predictions based on the smooth fitted to all sites
with(pred.aqs$pred, lines(date, smooth.fixed, col=1, lwd=2) )
with(pred.1site$pred, lines(date, smooth.fixed, col=2, lwd=2) )

##plot the predictions - One of the cases that vary in space, i.e. the smooth
##fit to the closest site.
##first extract as a data matrix
D <- with(pred.aqs$pred, createDataMatrix(obs=smooth.closest.fixed,
                                          date=date, ID=ID) )

##observations as a function of date
##(only five sites for clarity)
mesa.model <- dropObservations(mesa.model, !(mesa.model$obs$idx %in% c(1,2,3,23,24)))
plot(mesa.model, "loc.obs", type=as.factor(mesa.model$locations$ID),
     legend.loc=NULL, pch=19, cex=.25)
##Add the predictions based on the smooth
##fitted to the closest site
for(i in 1:5){
  lines(as.Date(rownames(D)), D[,mesa.model$locations$ID[i]], col=i, lwd=2)
}
