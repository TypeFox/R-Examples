## ---- echo=FALSE---------------------------------------------------------
library(credule)

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


## ---- echo=FALSE---------------------------------------------------------
t = c(1,3,5,7)
h = c(0.01,0.025,0.03,0.04)
spfun = function(i) {
  if (i==1) return(exp(-t[i]*h[i]))
  else return(spfun(i-1)*exp(-t[i]*h[i]))
}
sp = sapply(seq(1,4),spfun)
creditcurve_example = data.frame(tenor = as.character(t), hazardrate = percent(h), survprob = percent(sp))
knitr::kable(creditcurve_example, digits=4,caption = "Example of step-wise constant hazard rate function", row.names=FALSE)

## ---- fig.show='hold',echo=FALSE-----------------------------------------
barplot(h, names.arg = t, xlab="Time", ylab = "Non-Cummulative Hazard rate", main ="Hazard Rate (Non-Cum)")
plot(t, sp, type = "l", xlab="Time", ylab = "Survival Probability", main ="Survival Probability" )

## ---- echo=TRUE----------------------------------------------------------
yieldcurveTenor = c(1.0,2.0,3.0,4.0,5.0,7.0,
                    10.0,15.0,20.0,30.0)
yieldcurveRate = c(0.002585,0.005034,0.008981,0.012954,0.016452,
                   0.021811,0.027007,0.031718,0.033834,0.035056)

## ---- fig.show='hold',echo=FALSE-----------------------------------------
yieldcurve = data.frame(t(percent(yieldcurveRate)))
colnames(yieldcurve) = yieldcurveTenor
knitr::kable(yieldcurve,caption = "USD yield curve (%) as of 27 May 2014")

## ---- echo=TRUE----------------------------------------------------------
cdsTenor = c(1,2,3,4,5,7,10,15,20,30)
cdsSpread_PFE = c(0.0003,0.0009,0.0015,0.0021,0.0028,0.0043,0.0061,0.0063, 0.0068,0.0066)

## ---- fig.show='hold',echo=FALSE-----------------------------------------
cdsQuote_PFE = data.frame(t(cdsSpread_PFE*10000))
colnames(cdsQuote_PFE) = cdsTenor
knitr::kable(cdsQuote_PFE, digits=0,caption = "Pfizer (PFE) CDS spreads (basis points) as of 27 May 2014")

## ---- echo=TRUE----------------------------------------------------------
cdsTenor = c(1,2,3,4,5,7,10,15,20,30)
cdsSpread_RSH = c(0.6405,0.5956,0.5511,0.5144,0.4894,0.4511,0.4156,0.3815,0.3657,0.3506)

## ---- fig.show='hold',echo=FALSE-----------------------------------------
cdsQuote_RSH = data.frame(t(cdsSpread_RSH*10000))
colnames(cdsQuote_RSH) = cdsTenor
knitr::kable(cdsQuote_RSH, digits=0,caption = "Radioshack (RSH) CDS spreads (basis points) as of 27 May 2014")

## ---- echo=FALSE---------------------------------------------------------
# bootstrapping parameters
premiumFrequency = 4
defaultFrequency = 12
accruedPremium = TRUE
RR = 0.40

## ---- echo=TRUE----------------------------------------------------------
creditcurve_RSH = bootstrapCDS(yieldcurveTenor, yieldcurveRate, 
                   cdsTenor, cdsSpread_RSH,
                   RR, premiumFrequency, defaultFrequency, accruedPremium)

## ---- echo=FALSE---------------------------------------------------------
creditcurve_RSH_disp = creditcurve_RSH
creditcurve_RSH_disp$tenor = as.character(creditcurve_RSH_disp$tenor)
creditcurve_RSH_disp$hazrate = percent(creditcurve_RSH_disp$hazrate)
creditcurve_RSH_disp$survprob = percent(creditcurve_RSH_disp$survprob)
knitr::kable(creditcurve_RSH_disp, digits=4,caption = "RadioShack Corp Credit Curve", row.names=FALSE)

## ---- echo=FALSE, fig.height=3.5, fig.width=6----------------------------
default.par = par()
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
barplot(creditcurve_RSH$hazrate, names.arg = creditcurve_RSH$tenor, xlab="Time", ylab = "Non-Cummulative Hazard rate", main ="Hazard Rate (Non-Cum)")
plot(creditcurve_RSH$tenor, creditcurve_RSH$survprob, type = "l", xlab="Time", ylab = "Survival Probability", main ="Survival Probability",ylim=c(0,1))
mtext("RadioShack Corp on 27-May-2014", outer = TRUE, cex = 1.5)
par = default.par

## ---- echo=TRUE----------------------------------------------------------
creditcurve_PFE = bootstrapCDS(yieldcurveTenor, yieldcurveRate, 
                   cdsTenor, cdsSpread_PFE,
                   RR, premiumFrequency, defaultFrequency, accruedPremium)

## ---- echo=FALSE---------------------------------------------------------
creditcurve_PFE_disp = creditcurve_PFE
creditcurve_PFE_disp$tenor = as.character(creditcurve_PFE_disp$tenor)
creditcurve_PFE_disp$hazrate = percent(creditcurve_PFE_disp$hazrate)
creditcurve_PFE_disp$survprob = percent(creditcurve_PFE_disp$survprob)
knitr::kable(creditcurve_PFE_disp, digits=4,caption = "Pfizer Inc Credit Curve", row.names=FALSE)

## ---- echo=FALSE, fig.height=3.5, fig.width=6----------------------------
default.par = par()
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
barplot(creditcurve_PFE$hazrate, names.arg = creditcurve_PFE$tenor, xlab="Time", ylab = "Non-Cummulative Hazard rate", main ="Hazard Rate (Non-Cum)")
plot(creditcurve_PFE$tenor, creditcurve_PFE$survprob, type = "l", xlab="Time", ylab = "Survival Probability", main ="Survival Probability",ylim=c(0,1))
mtext("Pfizer Inc on 27-May-2014", outer = TRUE, cex = 1.5)
par = default.par

