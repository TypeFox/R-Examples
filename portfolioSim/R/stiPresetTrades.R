################################################################################
##
## $Id: stiPresetTrades.R 344 2006-10-01 05:06:05Z enos $
##
## Methods for class stiPresetTrades
##
################################################################################


setMethod("getSimTrades",
          signature(object = "stiPresetTrades",
                    period = "orderable",
                    holdings = "portfolio",
                    sim.data = "simData",
                    verbose = "logical"),
          function(object, period, holdings, sim.data, verbose){
            if(period %in% object@periods){                         
              trades <- object@sim.trades[object@periods == period][[1]]
            }
            else{
              trades <- new("trades")
            }

            return(new("simTrades", period = period, trades = trades))
          
          }
          )
