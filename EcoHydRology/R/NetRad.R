NetRad <- function (lat, Jday, Tx, Tn, albedo = 0.18, forest = 0, slope = 0, aspect = 0, airtemp = (Tn+Tx)/2, cloudiness = "Estimate", surfemissivity = 0.97, surftemp=(Tn+Tx)/2, units = "kJm2d", AEparams=list(vp=NULL, opt="linear")) {
                #  units : kJm3d or Wm2
        if (cloudiness[1] == "Estimate") {
                cloudiness <- EstCloudiness(Tx, Tn)
        }
    if (units == "kJm2d") convert <- 1 else convert <- 86.4
        return( signif((Solar(lat, Jday, Tx, Tn, albedo, forest, slope, aspect, printWarn=FALSE) + Longwave(AtmosphericEmissivity(airtemp, cloudiness, AEparams$vp, AEparams$opt),airtemp) - Longwave(surfemissivity, surftemp)) / convert , 5) )
}
