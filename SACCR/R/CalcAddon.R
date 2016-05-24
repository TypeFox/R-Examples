#' Calculates the aggregate amount of the addon after splitting per asset class
#' and dividing the trades into the corresponding netting sets per currency, timebucket  etc.
#' 
#' @title Calculates the Addon amount
#' @param trades The full list of the Trade Objects
#' @param MF (Optional) The Maturity Factor based on the collateral agreement  
#' @param factor_mult (Optional) The Multiplication Factor applicable for volatility/basis trades  
#' @return The aggregate amount of the addon summed up for all the asset classes
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
#' 
CalcAddon <- function(trades, MF, factor_mult)  {
## function which will calculate the Add-On for all the trade classes
  
  superv <- LoadSupervisoryData()
  
  if(missing(factor_mult))
    factor_mult = 1
  
  # finding all the different trade classes that have been sent to the function
  trade_classes <- unique(lapply(trades, function(x) x$TradeGroup))
  # this array will hold the addon for all the trade classes
  trade_classes_addon <- array(data<-0,dim<-length(trade_classes))
  
  # going through each trade class
  for (i in 1:length(trade_classes))
  {      
    #picking up the trades belonging to this specific trade class
    group_trades <- trades[sapply(trades, function(x) x$TradeGroup==trade_classes[i])]
    
    
    if(trade_classes[i]=='FX')
    {
      # for the FX case the Hedging sets will be created based on the ccy pair
      ccypairs   <- unique(lapply(group_trades, function(x) x$ccyPair))
      ccypairs_addon <- array(data<-0,dim<-length(ccypairs))
      
      # going through all the ccy pairs found
      for (j  in 1:length(ccypairs))
      {  
        # picking up all the trades for this ccy pair
        ccypairs_trades  <- group_trades[sapply(group_trades, function(x) (x$ccyPair==ccypairs[j]))]
        
        # for each of the trade calculate the Adjusted Notional and their contribution to the addon of the hedging set
        for (l in 1:length(ccypairs_trades))
        {
          AdjNotional         <- ccypairs_trades[[l]]$CalcAdjNotional()
          
          if(missing(MF))
          {
          # calculate maturity factor
          maturity_factor <- ccypairs_trades[[l]]$CalcMaturityFactor()
          } else
          { maturity_factor = MF}
          # if the trade is option based then for the delta calculation the volatility will be used
          if (ccypairs_trades[[l]]$TradeType=='Option')
          {
            volatility   <- superv$Supervisory_option_volatility[superv$Asset_Class==ccypairs_trades[[l]]$TradeGroup&superv$SubClass==ccypairs_trades[[l]]$SubClass]
            superv_delta <- ccypairs_trades[[l]]$CalcSupervDelta(volatility)
          }
          else
          {
            superv_delta <- ccypairs_trades[[l]]$CalcSupervDelta()
          }
          # aggregating the add-on contribution for a specific hedging set
          ccypairs_addon[j] <- ccypairs_addon[j] + superv_delta*AdjNotional*maturity_factor
        }
      }
      # getting the supervisory factor
      supervisory_factor <- factor_mult*superv$Supervisory_factor[superv$Asset_Class==ccypairs_trades[[l]]$TradeGroup&superv$SubClass==ccypairs_trades[[l]]$SubClass]
      
      # adding up the addon of the hedging set to the trade class
      trade_classes_addon[i]<- trade_classes_addon[i] + supervisory_factor*sum(abs(ccypairs_addon))
    }  else if(trade_classes[i]=='IRD')
    {
      # setting the time bucket for each of the trades
      lapply(group_trades, function(x) x$TimeBucket <- x$SetTimeBucket())
      # picking up the currencies found in the IRD trades which will be the first-level grouping applied 
      currencies   <- unique(lapply(group_trades, function(x) x$Currency))
      currencies_addon <- array(data<-0,dim<-length(currencies))
      for (j  in 1:length(currencies))
      {  
        currency_trades  <- group_trades[sapply(group_trades, function(x) (x$Currency==currencies[j]))]
        
        # after picking up the trades related to a specific currency, a second-level grouping will take place
        # based on the time buckets
        timebuckets       <- unique(lapply(currency_trades, function(x) x$TimeBucket))
        timebuckets_addon <- array(data<-0,dim<-3)
        
        for (k in 1:length(timebuckets))
        {  
          #picking up all the trades belonging to a specific timebucket
          timebuckets_trades  <- currency_trades[sapply(currency_trades, function(x) x$TimeBucket==timebuckets[k])]
          
          for (l in 1:length(timebuckets_trades))
          {
            if(missing(MF))
            {
            # calculate maturity factor
            maturity_factor     <- timebuckets_trades[[l]]$CalcMaturityFactor()
            } else
            {
              maturity_factor = MF
            }
            AdjNotional         <- timebuckets_trades[[l]]$CalcAdjNotional()
            if (timebuckets_trades[[l]]$TradeType=='Option')
            {
              volatility <- superv$Supervisory_option_volatility[superv$Asset_Class==timebuckets_trades[[l]]$TradeGroup&superv$SubClass==timebuckets_trades[[l]]$SubClass]
              superv_delta        <- timebuckets_trades[[l]]$CalcSupervDelta(volatility)
            }
            else
            {
              superv_delta <- timebuckets_trades[[l]]$CalcSupervDelta()
            }
            timebuckets_addon[timebuckets[[k]]] <- timebuckets_addon[timebuckets[[k]]] + superv_delta*AdjNotional*maturity_factor
          }
        }
        # aggregating the add-on timebuckets recognizing correlation between each time bucket  
        currencies_addon[j] <- (timebuckets_addon[1]^2+timebuckets_addon[2]^2+timebuckets_addon[3]^2+1.4*timebuckets_addon[2]*timebuckets_addon[3]+1.4*timebuckets_addon[2]*timebuckets_addon[1]+0.6*timebuckets_addon[2]*timebuckets_addon[1])^0.5
        supervisory_factor <- factor_mult*superv$Supervisory_factor[superv$Asset_Class==timebuckets_trades[[l]]$TradeGroup&superv$SubClass==timebuckets_trades[[l]]$SubClass]
        
        # adding up the addon of each currency after multiplying with the supervisory factor
        trade_classes_addon[i]<- trade_classes_addon[i] + supervisory_factor*currencies_addon[j]
      }
    }else  if(trade_classes[i]=='Credit')
    {
      # for the Credit case the Hedging sets will be created based on the reference entity
      refEntities   <- unique(lapply(group_trades, function(x) x$RefEntity))
      refEntities_addon <- array(data<-0,dim<-length(refEntities))
      supervisory_corel <- array(data<-0,dim<-length(refEntities))
      for (j  in 1:length(refEntities))
      {  
        refEntities_trades  <- group_trades[sapply(group_trades, function(x) (x$RefEntity == refEntities[j]))]
        
        for (k in 1:length(refEntities_trades))
        { 
          if(missing(MF))
          {
          # calculate maturity factor
          maturity_factor     <- refEntities_trades[[k]]$CalcMaturityFactor()
          } else
          { maturity_factor <- MF }
          AdjNotional         <- refEntities_trades[[k]]$CalcAdjNotional()
          superv_delta        <- refEntities_trades[[k]]$CalcSupervDelta()
          refEntities_addon[j]<- refEntities_addon[j] + AdjNotional*superv_delta*maturity_factor
        }
        
        AssetClass<-paste(refEntities_trades[[1]]$TradeGroup,refEntities_trades[[1]]$TradeType,sep="")
        supervisory_factor <- factor_mult*superv$Supervisory_factor[superv$Asset_Class==AssetClass&superv$SubClass==refEntities_trades[[1]]$SubClass]
        refEntities_addon[j] <- refEntities_addon[j]*supervisory_factor
        supervisory_corel[j]  <- superv$Correlation[superv$Asset_Class==AssetClass&superv$SubClass==refEntities_trades[[1]]$SubClass]
      }
      systematic_component     <- (sum(refEntities_addon*supervisory_corel))^2
      idiosynchratic_component <-  sum((rep(1,length(refEntities))-supervisory_corel^2)*refEntities_addon^2)
      trade_classes_addon[i]   <- sqrt(systematic_component + idiosynchratic_component)
    }else   if(trade_classes[i]=='Commodity')
    {
      AssetClass <-'Commodity'
      
      # for the commodities case the Hedging sets will be created based on the commodity sector on a first
      # level (energy,metals etc) and on a second level based on the actual commodities (oil, silver etc)
      HedgingSets   <- unique(lapply(group_trades, function(x) x$SubClass))
      HedgingSets_addon <- array(data<-0,dim<-length(HedgingSets))
      
      for (j  in 1:length(HedgingSets))
      {  
        HedgingSets_trades  <- group_trades[sapply(group_trades, function(x) (x$SubClass == HedgingSets[j]))]
        
        com_types       <- unique(lapply(HedgingSets_trades, function(x) x$commodity_type))
        com_types_addon <- array(data<-0,dim<-length(com_types))
        
        for (k in 1:length(com_types))
        {  
          com_types_trades  <- HedgingSets_trades[sapply(HedgingSets_trades, function(x) x$commodity_type==com_types[k])]
          
          
          for (l in 1:length(com_types_trades))
          { 
            if(missing(MF))
            {
            # calculate maturity factor
            maturity_factor     <- com_types_trades[[l]]$CalcMaturityFactor()
            } else
            { maturity_factor = MF}
            AdjNotional         <- com_types_trades[[l]]$CalcAdjNotional()
            superv_delta        <- com_types_trades[[l]]$CalcSupervDelta()
            com_types_addon[k]<- com_types_addon[k] + AdjNotional*superv_delta*maturity_factor
          }
          
          supervisory_factor <- factor_mult*superv$Supervisory_factor[superv$Asset_Class==AssetClass&(superv$SubClass==HedgingSets_trades[[1]]$SubClass|superv$SubClass==HedgingSets_trades[[1]]$commodity_type)]
          com_types_addon[k]<- com_types_addon[k]*supervisory_factor
        }
        supervisory_corel     <- superv$Correlation[superv$Asset_Class==AssetClass&(superv$SubClass==HedgingSets_trades[[1]]$SubClass|superv$SubClass==HedgingSets_trades[[1]]$commodity_type)]
        HedgingSets_addon[j]  <- sqrt((sum(com_types_addon)*supervisory_corel)^2 + (1-supervisory_corel^2)*sum(com_types_addon^2))
        
      }
      trade_classes_addon[i]   <- sum(HedgingSets_addon)
    }  
  }
  return(sum(trade_classes_addon))

}