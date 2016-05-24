obese <-
function(prevalence,month,weight_gain,path = system.file("OB_above.RData", "OB_below.RData", package = "GWG"))
{
data(OB_above)
data(OB_below)
i <- month
PPVex  <- (OB_above$se.est[OB_above$month==i]*prevalence)/ (OB_above$se.est[OB_above$month==i]*prevalence+(1-OB_above$sp.est[OB_above$month==i])*(1-prevalence))
PPVin  <- (OB_below$se.est[OB_below$month==i]*prevalence)/ (OB_below$se.est[OB_above$month==i]*prevalence+(1-OB_below$sp.est[OB_below$month==i])*(1-prevalence))
NPVex   <-(OB_above$sp.est[OB_above$month==i]*(1-prevalence))/ (OB_above$sp.est[OB_above$month==i]*(1-prevalence)+(1-OB_above$se.est[OB_above$month==i])*(prevalence))
NPVin  <- (OB_below$sp.est[OB_below$month==i]*(1-prevalence))/ (OB_below$sp.est[OB_above$month==i]*(1-prevalence)+(1-OB_below$se.est[OB_below$month==i])*(prevalence))
valuesex <- list(month=month, prevalence=prevalence, PPVexcessive=PPVex,NPVexcessive=NPVex)
valuesin <- list(month=month, prevalence=prevalence, PPVinadequate=PPVin,NPVinadequate=NPVin)
if (weight_gain == "excessive") {
        return(valuesex)
    }
if (weight_gain == "inadequate") {
        return(valuesin)
    }
}

