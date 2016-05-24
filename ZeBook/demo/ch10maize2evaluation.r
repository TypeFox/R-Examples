################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-14
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
list_n_sy=unique(maize.data_EuropeEU$sy)

# 1) simulations for our n sites
sim = maize.multisy(maize.define.param()["nominal",], list_n_sy, 100,250, weather_all=weather_EuropeEU)
# 2) Evaluation - combining simulation and data in a single table
# table of simulation corresponding to available observations
simobs <- merge(sim,maize.data_EuropeEU, by=c("sy",  "day"), all.x=TRUE)
# table of simulation corresponding to observations of biomass at day=240
simobs240<-subset(simobs, day==240)
simobs240

# 3) Evaluation - plot simulation and observation on a same graph
par(mfrow=c(5,3), mar=c(2,2,1,0.5))
for(sy in list_n_sy){
  result <- simobs[simobs$sy==sy,]
  plot(result$day,result$B,type="l",xlim=c(100,250),ylim=c(0,4000))
  title (sy)
  points(result$day,result$Bobs)
}

# 4) Evaluation - plot observation versus simulation for biomass at day=240
dev.new()
par(mfrow=c(1,1))
plot(simobs240$B,simobs240$Bobs, xlim=range(c(simobs240$B,simobs240$Bobs)),ylim=range(c(simobs240$B,simobs240$Bobs)), ylab="Bobs", xlab="B")
abline(a=0,b=1)
dev.new()
plot(simobs240$B,simobs240$Bobs-simobs240$B, xlab="predicted values",ylab="residuals")
abline(a=0,b=0)

# 5) Evaluation - criteria of goodness of fit
goodness.of.fit(simobs240$B,simobs240$Bobs,draw.plot=FALSE)

# end of file