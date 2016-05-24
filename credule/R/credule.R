priceCDS = function(yieldcurveTenor,
                                             yieldcurveRate,
                                             creditcurveTenor,
                                             creditcurveSP,
                                             cdsTenors,
                                             recoveryRate,
                                             numberPremiumPerYear = c(4,2,1,12),
                                             numberDefaultIntervalPerYear = 12,
                                             accruedPremium = c(TRUE,FALSE)) 
  {
  if (length(yieldcurveTenor) != length(yieldcurveRate)) {
    warning("yieldcurveTenor and yieldcurveRate do not have the same length")
  }
  else if (length(creditcurveTenor) != length(creditcurveSP)) {
    warning("creditcurveTenor and creditcurveSP do not have the same length")
  }
  else {
    result = .C("priceCreditDefaultSwapSpreads", 
                yieldcurve = as.numeric(c(yieldcurveTenor,yieldcurveRate)), 
                nyieldcurve = as.integer(length(yieldcurveTenor)),
                creditcurve = as.numeric(c(creditcurveTenor,creditcurveSP)),
                ncreditcurve = as.integer(length(creditcurveTenor)),
                cdsTenors = as.numeric(cdsTenors),
                ncdsTenors = as.integer(length(cdsTenors)),
                numberPremiumPerYear = as.integer(numberPremiumPerYear),
                numberDefaultIntervalPerYear = as.integer(numberDefaultIntervalPerYear),
                accruedPremiumFlag = as.integer(accruedPremium),
                recoveryRate = as.numeric(recoveryRate),
                spreads = as.numeric(rep(0,length(cdsTenors))),
                warningFlag = as.numeric(0)
    )
    
    df = data.frame(tenor=cdsTenors,spread=result$spreads)
    
    return (df)               
  }  
}

bootstrapCDS = function(yieldcurveTenor,
                        yieldcurveRate,
                        cdsTenors,
                        cdsSpreads,
                        recoveryRate,
                        numberPremiumPerYear = c(4,2,1,12),
                        numberDefaultIntervalPerYear = 12,
                        accruedPremium = c(TRUE,FALSE)) 
  {
  if (length(yieldcurveTenor) != length(yieldcurveRate)) {
    warning("yieldcurveTenor and yieldcurveRate do not have the same length")
  }
  else if (length(cdsTenors) != length(cdsSpreads)) {
    warning("cdsTenors and cdsSpreads do not have the same length")
  }
  else {
    result = .C("bootstrapCreditDefaultSwapSpreads", 
                yieldcurve = as.numeric(c(yieldcurveTenor,yieldcurveRate)), 
                nyieldcurve = as.integer(length(yieldcurveTenor)),
                cdsTenors = as.numeric(cdsTenors),
                ncdsTenors = as.integer(length(cdsTenors)),
                spreads = as.numeric(cdsSpreads),
                numberPremiumPerYear = as.integer(numberPremiumPerYear),
                numberDefaultIntervalPerYear = as.integer(numberDefaultIntervalPerYear),
                accruedPremiumFlag = as.integer(accruedPremium),
                recoveryRate = as.numeric(recoveryRate),
                output = as.numeric(rep(0,length(cdsTenors)*2)),
                warningFlag = as.numeric(0)
    )
    
    if (result$warningFlag) {
      stop( "CDS term structure can't be bootstrapped. Input CDS spreads would imply increasing survival probabilities over time" )
    }
    
    df = data.frame(tenor=cdsTenors,survprob=result$output[1:length(cdsTenors)],hazrate=result$output[(length(cdsTenors)+1):(2*length(cdsTenors))])
    
    
    return (df)               
  }  
}