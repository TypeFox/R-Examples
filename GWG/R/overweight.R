overweight <-
function(prevalence,month,weight_gain,path = system.file("OW_above.RData", "OW_below.RData", package = "GWG"))
{
data(OW_above)
data(OW_below)
i <- month
PPVex  <- (OW_above$se.est[OW_above$month==i]*prevalence)/ (OW_above$se.est[OW_above$month==i]*prevalence+(1-OW_above$sp.est[OW_above$month==i])*(1-prevalence))
PPVin  <- (OW_below$se.est[OW_below$month==i]*prevalence)/ (OW_below$se.est[OW_above$month==i]*prevalence+(1-OW_below$sp.est[OW_below$month==i])*(1-prevalence))
NPVex   <-(OW_above$sp.est[OW_above$month==i]*(1-prevalence))/ (OW_above$sp.est[OW_above$month==i]*(1-prevalence)+(1-OW_above$se.est[OW_above$month==i])*(prevalence))
NPVin  <- (OW_below$sp.est[OW_below$month==i]*(1-prevalence))/ (OW_below$sp.est[OW_above$month==i]*(1-prevalence)+(1-OW_below$se.est[OW_below$month==i])*(prevalence))
valuesex <- list(month=month, prevalence=prevalence, PPVexcessive=PPVex,NPVexcessive=NPVex)
valuesin <- list(month=month, prevalence=prevalence, PPVinadequate=PPVin,NPVinadequate=NPVin)
if (weight_gain == "excessive") {
        return(valuesex)
    }
if (weight_gain == "inadequate") {
        return(valuesin)
    }
}

