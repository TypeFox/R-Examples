############################################


setGeneric(
  name="Asset_PriceDistribution",
  def = function(.Object, type, t,T,nCoupons=NULL,couponsRate=NULL,omega=NULL, s=NULL, Strike=NULL)
  {
    standardGeneric("Asset_PriceDistribution")
  }
)

setMethod(
  f="Asset_PriceDistribution",
  signature="Scenarios",
  definition=function(.Object, type, t, T, nCoupons=NULL, couponsRate=NULL, omega=NULL, s=NULL, Strike=NULL)
  {
    if (missing(type)||!(type %in% c("Zero-Coupon", "Bond", "CBond", "ConvBond", "EuroCall_Stock", "EuroPut_Stock", "EuroCall_ZC", "EuroPut_ZC", "CDS"))) 
    {stop("Type argument must be : 'Zero-Coupon', 'Bond', 'CBond', 'ConvBond', 'EuroCall_Stock', 'EuroPut_Stock', 'EuroCall_ZC', 'EuroPut_ZC', 'CDS' ")}
    
    if (missing(t)||missing(T)) {stop("You must fill at least the parameters t and T.")} else
    {
      
      y <- 0
      
      if (type == "Zero-Coupon")  
      {
        if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)||!missing(s)||!missing(Strike)) warning("Only these parameters are actually used for this asset : t, T")
        y <- ZCBond_PriceDistribution(.Object,t,T)
      }
      
      if (type == "Bond") 
      {    
        if (missing(nCoupons)||missing(couponsRate)) {stop("You must fill the parameters nCoupons and couponsRate")} else {
          if (!missing(omega)||!missing(s)||!missing(Strike)) warning("Only these parameters are actually used for this asset : t, T, nCoupons, couponsRate")
          y <- Bond_PriceDistribution(.Object,t,T,nCoupons,couponsRate)}
      }
      
      if (type == "CBond") 
      {
        if (missing(nCoupons)||missing(couponsRate)||missing(omega)) {stop("You must fill the parameters nCoupons, couponsRate and omega.")} else {
          if (!missing(s)||!missing(Strike)) warning("Only these parameters are actually used for this asset : t,T,nCoupons,couponsRate,omega")
          y <- CBond_PriceDistribution(.Object,t,T,nCoupons,couponsRate,omega)}
      }
      
      if (type == "EuroCall_Stock")  
      {
        if (missing(Strike)) {stop("You must fill the Strike parameter")} else {
          if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)||!missing(s)) warning("Only these parameters are actually used for this asset: t, T, Strike")
          y <- EuroCall_Stock_PriceDistribution(.Object,t,T,Strike)}
      }
      
      if (type == "EuroPut_Stock")  
      {
        if (missing(Strike)) {stop("You must fill the Strike parameter")} else {	
          if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)||!missing(s)) warning("Only these parameters are actually used for this asset : t, T, Strike")
          y <- EuroPut_Stock_PriceDistribution(.Object,t,T,Strike)}
      }
      
      if (type == "EuroCall_ZC")  
      { 
        if (missing(s)||missing(Strike)) {stop("You must fill the s and Strike parameters")} else {
          if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)) warning("Only these parameters are actually used for this asset : t, T, s, Strike")
          y <- EuroCall_ZC_PriceDistribution(.Object,t,T,s,Strike)}
      }
      
      if (type == "EuroPut_ZC")  
      {
        if (missing(s)||missing(Strike)) {stop("You must fill the s and Strike parameters")} else {
          if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)) warning("Only these parameters are actually used for this asset : t, T, s, Strike") 
          y <- EuroPut_ZC_PriceDistribution(.Object,t,T,s,Strike)}
      }
      
      if (type == "CDS")  
      {
        if (missing(omega)||missing(t)||missing(T)) {stop("You must fill the s and Strike parameters")} else {
          if (!missing(nCoupons)||!missing(couponsRate)||!missing(Strike)) warning("Only these parameters are actually used for this asset : omega, t, T") 
          y <- CDSPremium_PriceDistribution(.Object,omega=omega,t=t,T=T)}
      }
      
      if (type == "ConvBond") 
      {    
        if (missing(nCoupons)||missing(couponsRate)) {stop("You must fill the parameters nCoupons and couponsRate")} else {
          if (!missing(omega)||!missing(s)||!missing(Strike)) warning("Only these parameters are actually used for this asset : t, T, nCoupons, couponsRate")
          y <- ConvBond_PriceDistribution(.Object,t,T,nCoupons,couponsRate)}
      }
      
      return (y)
      
    }
  }
)