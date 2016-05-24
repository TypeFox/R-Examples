### R code from vignette source 'ssfavignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: install
###################################################
#install.packages("ssfa") 
library(ssfa) 


###################################################
### code chunk number 2: data
###################################################
data(SSFA_example_data)
data(Italian_W)
names(SSFA_example_data)


###################################################
### code chunk number 3: fig1
###################################################
plot(SSFA_example_data$log_x, SSFA_example_data$log_y, pch=16, cex=0.5, xlab="log_x", ylab="log_y", cex.axis=0.8, main="SSFA_example_data")


###################################################
### code chunk number 4: sfa
###################################################
sfa <- ssfa(log_y ~ log_x , data = SSFA_example_data, data_w=Italian_W, 
            form = "production", par_rho=FALSE)
summary(sfa)


###################################################
### code chunk number 5: moran_code1
###################################################
 moran.test(residuals(sfa), listw=sfa$list_w)


###################################################
### code chunk number 6: moran_code2
###################################################
 plot_moran(sfa, listw=sfa$list_w)


###################################################
### code chunk number 7: fig2
###################################################
 plot_moran(sfa, listw=sfa$list_w)


###################################################
### code chunk number 8: ssfa
###################################################
ssfa <- ssfa(log_y ~ log_x , data = SSFA_example_data, data_w=Italian_W, 
             form = "production", par_rho=TRUE)
summary(ssfa)


###################################################
### code chunk number 9: moran_code3
###################################################
 moran.test(residuals(ssfa), listw=ssfa$list_w)


###################################################
### code chunk number 10: moran_code4
###################################################
 plot_moran(ssfa, listw=sfa$list_w)


###################################################
### code chunk number 11: fig3
###################################################
 plot_moran(ssfa, listw=sfa$list_w)


###################################################
### code chunk number 12: fitted
###################################################
ssfa_fitted <- fitted.ssfa(ssfa)
sfa_fitted <- fitted.ssfa(sfa)


###################################################
### code chunk number 13: plot_code
###################################################
plot_fitted(SSFA_example_data$log_x, SSFA_example_data$log_y, ssfa, pch=16, cex=0.5,
              xlab="X", ylab="Y", cex.axis=0.8 )
points(SSFA_example_data$log_x, SSFA_example_data$log_y, pch=16, cex=0.5, 
  col= ifelse(eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 0.20) , "#D7191C",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.20)
          &eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 0.4) ,"#FF8C00",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.4)
          &eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 0.6) ,"#FFFF00",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.6) 
          &eff.ssfa(ssfa)<quantile(eff.ssfa(ssfa), 0.8) ,"#ADFF2F",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.8) 
          &eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 1),"#008B00", "#2F4F4F"))))))
lines(sort(SSFA_example_data$log_x),sfa_fitted[order(SSFA_example_data$log_x)],
            col="red")


###################################################
### code chunk number 14: fig4
###################################################
plot_fitted(SSFA_example_data$log_x, SSFA_example_data$log_y, ssfa, pch=16, cex=0.5, xlab="X", ylab="Y", cex.axis=0.8 )
points(SSFA_example_data$log_x, SSFA_example_data$log_y, pch=16, cex=0.5, col=ifelse(eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 0.20) , "#D7191C",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.20)&eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 0.4) ,"#FF8C00",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.4) &eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 0.6) ,"#FFFF00",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.6) &eff.ssfa(ssfa)<quantile(eff.ssfa(ssfa), 0.8) ,"#ADFF2F",
        ifelse(eff.ssfa(ssfa)>quantile(eff.ssfa(ssfa), 0.8) &eff.ssfa(ssfa)<=quantile(eff.ssfa(ssfa), 1),"#008B00","#2F4F4F"))))))
lines(sort(SSFA_example_data$log_x), sfa_fitted[order(SSFA_example_data$log_x)],col="red")


###################################################
### code chunk number 15: residuals
###################################################
ssfa_residuals <- residuals.ssfa(ssfa)
sfa_residuals <- residuals.ssfa(sfa)


###################################################
### code chunk number 16: efficiency
###################################################
ssfa_eff <-  eff.ssfa(ssfa)
#sfa_eff <-  eff.ssfa(sfa)

#summary(sfa_eff)
#summary(ssfa_eff)

ssfa_u <- u.ssfa(ssfa)
#sfa_u <- u.ssfa(sfa)

#summary(ssfa_u)
#summary(sfa_u)


