# The Observed Properties of Liquid Helium at the Saturated Vapor Pressure

# Russell J. Donnelly, 2004
# The Observed Properties of Liquid Helium at the Saturated Vapor Pressure
# http://darkwing.uoregon.edu/~rjd/vapor17.htm
# Figure 16.1. The recommended values for the latent heat of vaporization of liquid 4He as a function of temperature at the saturated vapor pressure.
plot(recommendedLatentHeatOfVaporizationOfLiquidHe4[,1], recommendedLatentHeatOfVaporizationOfLiquidHe4[,2],type='l',xlab='Temperature(K)',ylab='L(J/mol)')
# Figure 16.2. Detail of the recommended values latent heat of vaporization about the lambda transition.
plot(log(abs(recommendedLatentHeatOfVaporizationOfLiquidHe4[,1]-2.1768)), recommendedLatentHeatOfVaporizationOfLiquidHe4[,2],type='l',xlab='Log|T-T{lambda}|',ylab='L(J/mol)')

# Figure 16.3. The fractional deviation of values of the adopted database from the recommended values for the latent heat expressed in percent.
f <- approxfun(recommendedLatentHeatOfVaporizationOfLiquidHe4[,1], recommendedLatentHeatOfVaporizationOfLiquidHe4[,2])
r <- f(adoptedLatentHeatOfVaporizationOfLiquidHe4[,1])/adoptedLatentHeatOfVaporizationOfLiquidHe4[,2]
plot(adoptedLatentHeatOfVaporizationOfLiquidHe4[,1], r,type='l')