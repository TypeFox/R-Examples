######################################################################
# Create the base Vol class
#
# This is used to represent the swap product and it will contained to all the swap-like classes

Vol = setRefClass("Vol",
                  fields = list(vol_strike   = "numeric",
                                annualization_factor = "numeric",
                                vega_notional = "numeric",
                                reference = "character"
                  ),
                  methods = list(
                    ComputeVarianceUnits = function() {
                     return(vega_notional/(2*vol_strike))
                    }
                  ))