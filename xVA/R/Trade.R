######################################################################
# Create the base Trade class
#
# This class contains fields/methods that are used from all the derived trade objects
Trade = setRefClass("Trade",
                     fields = list(     Notional   = "numeric",
                                        MtM        = "numeric",
                                        Currency   = "character",
                                        Si         = "numeric",
                                        Ei         = "numeric",
                                        BuySell    = "character",
                                        TradeGroup = "character",
                                        TradeType  = "character",
                                        SubClass   = "character"
                                        ),
                     methods = list(
                       CalcAdjNotional = function() {
                         ## calculates the adjusted notional by multiplying the notional amount with the
                         ## supervisory duration
                         if (TradeGroup=='IRD'||TradeGroup=='Credit')
                         {
                           AdjustedNotional = Notional * CalcSupervDuration() ;
                         } else
                         {
                           AdjustedNotional = Notional;
                         }
                         
                         return(AdjustedNotional)
                       },
                       #' @description Test Test
                       CalcSupervDuration = function() {
                         ## calculates the supervisory duration (applicable for IRDs and Credit derivatives)
                         if(Ei<1)
                         {SupervisoryDuration = sqrt(Ei)
                         }else SupervisoryDuration = (exp(-0.05*Si)-exp(-0.05*Ei))/0.05;
                         
                         return(SupervisoryDuration)
                       },
                       CalcMaturityFactor = function() {
                         ## calculates the maturity factor 
                         if(Ei<1)
                         {MaturityFactor = sqrt(Ei)
                         }else MaturityFactor = 1;
                         
                         return(MaturityFactor)
                       },
                       CalcSupervDelta = function(Superv_Vol) {
                         if (missing(Superv_Vol))
                         {
                           if(BuySell=="Buy")   return(1)
                           if(BuySell=="Sell")  return(-1)
                         } else
                         {
                           # in the option case the supervisory volatility is being populated
                           # the delta calculation is based on the Black-Scholes formula  
                           temp = (log(UnderlyingPrice/StrikePrice)+0.5*Superv_Vol^2*Si)/Superv_Vol*Si^0.5;
                           
                           if(BuySell=="Buy")
                           {
                             if(OptionType=='Call') return(pnorm(temp))
                             if(OptionType=='Put')  return(-pnorm(-temp))
                           }
                           if(BuySell=="Sell") 
                           {
                             if(OptionType=='Call') return(-pnorm(temp))
                             if(OptionType=='Put')  return(pnorm(-temp))
                           }
                         }
                       }
                     ))