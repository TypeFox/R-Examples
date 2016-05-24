##cross-validation load data
data(est.cv.mesa)
##...and old estimates
data(est.mesa.model)
##estimated parameters
par.cov <- coef(est.mesa.model, "cov")
par.all <- coef(est.mesa.model)
 
##boxplot of the different estimates from the CV
par(mfrow=c(1,1), mar=c(13,2.5,2,.5), las=2)
boxplot( est.cv.mesa, plot.type="cov", boxwex=.5)
##compare with estimates for all data
points((1:length(par.cov$par))+.35, par.cov$par, pch=19, col=2)
##and uncertainties
for(i in 1:length(par.cov$par)){
  lines(c(i,i)+.35, par.cov$par[i]+c(-1,1)*1.96*par.cov$sd[i], col=2, lwd=2)
}

##For all the parameters but with offset lines
par(mfrow=c(1,1), mar=c(13,2.5,2,.5), las=2)
boxplot(est.cv.mesa, plot.type="all", boxwex=.4, col="grey",
        main="Cross-validation estimates")
##compare with estimates for all data
points((1:length(par.all$par))+.35, par.all$par, pch=19, col=2)
##and uncertainties
for(i in 1:length(par.all$par)){
  lines(c(i,i)+.35, par.all$par[i]+c(-1,1)*1.96*par.all$sd[i], col=2, lwd=2)
}
