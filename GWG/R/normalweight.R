normalweight <-
function(prevalence,month,weight_gain,path = system.file("NW_above.RData", "NW_below.RData", package = "GWG"))
{
i <- month
PPVex  <- (NW_above$se.est[NW_above$month==i]*prevalence)/ (NW_above$se.est[NW_above$month==i]*prevalence+(1-NW_above$sp.est[NW_above$month==i])*(1-prevalence))
PPVin  <- (NW_below$se.est[NW_below$month==i]*prevalence)/ (NW_below$se.est[NW_above$month==i]*prevalence+(1-NW_below$sp.est[NW_below$month==i])*(1-prevalence))
NPVex   <-(NW_above$sp.est[NW_above$month==i]*(1-prevalence))/ (NW_above$sp.est[NW_above$month==i]*(1-prevalence)+(1-NW_above$se.est[NW_above$month==i])*(prevalence))
NPVin  <- (NW_below$sp.est[NW_below$month==i]*(1-prevalence))/ (NW_below$sp.est[NW_above$month==i]*(1-prevalence)+(1-NW_below$se.est[NW_below$month==i])*(prevalence))
valuesex <- list(month=month, prevalence=prevalence, PPVexcessive=PPVex,NPVexcessive=NPVex)
valuesin <- list(month=month, prevalence=prevalence, PPVinadequate=PPVin,NPVinadequate=NPVin)
if (weight_gain == "excessive") {
        return(valuesex)
    }
if (weight_gain == "inadequate") {
        return(valuesin)
    }
}
