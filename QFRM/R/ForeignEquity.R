#' @title ForeignEquity option valuation via Black-Scholes (BS) model
#' @description ForeignEquity Option via Black-Scholes (BS) model
#' @author Chengwei Ge, Department of Statistics, Rice University, 2015
#' 
#' @param o An object of class \code{OptPx}
#' @param Type ForeignEquity option type: 'Foreign' or 'Domestic'
#' @param I1 A spot price of the underlying security 1 (usually I1)
#' @param I2 A spot price of the underlying security 2 (usually I2)
#' @param g1 is the payout rate of the first stock
#' @param sigma1 a vector of implied volatilities for the associated security 1
#' @param sigma2 a vector of implied volatilities for the associated security 2 
#' @param rho is the correlation between asset 1 and asset 2
#' @details 
#' Two types of ForeignEquity options are priced: \code{'Foreign'} and \code{'Domestic'}. 
#' See "Exotic Options", 2nd, Peter G. Zhang for more details. 
#'  
#' @return A list of class \code{ForeignEquityBS} consisting of the original \code{OptPx} object 
#' and the option pricing parameters \code{I1},\code{I2}, \code{Type}, \code{isForeign}, and \code{isDomestic}
#' as well as the computed price \code{PxBS}.
#' @references Zhang, Peter G.  \emph{Exotic Options}, 2nd, 1998.  
#' 
#' @examples
#' o = OptPx(Opt(Style = 'ForeignEquity', Right = "Put"), r= 0.03)
#' ForeignEquityBS(o, I1=1540, I2=1/90, g1=.02, sigma1=.14,sigma2=0.18, rho=.03,Type='Foreign')
#' 
#' o = OptPx(Opt(Style = 'ForeignEquity',  Right = "Put", ttm=9/12, K=1600), r=.03)
#' ForeignEquityBS(o, I1=1540, I2=1/90, g1=.02, sigma1=.14,sigma2=0.18, rho=0.03,Type='Foreign')
#' 
#' o = OptPx(Opt(Style = 'ForeignEquity', Right = "C", ttm=9/12, K=1600), r=.03)
#' ForeignEquityBS(o, I1=1540, I2=1/90, g1=.02, sigma1=.14,sigma2=0.18, rho=0.03,Type='Foreign')
#' 
#' o = OptPx(Opt(Style = 'ForeignEquity', Right = "C", ttm=9/12, K=1600), r=.03)
#' ForeignEquityBS(o, I1=1540, I2=1/90, g1=.02, sigma1=.14,sigma2=0.18, rho=0.03,Type='Domestic')
#' 
#' o = OptPx(Opt(Style = 'ForeignEquity', Right = "P", ttm=9/12, K=1600), r=.03)
#' ForeignEquityBS(o, I1=1540, I2=1/90, g1=.02, sigma1=.14,sigma2=0.18, rho=0.03,Type='Domestic')
#' @export
#' 
ForeignEquityBS=function(o=OptPx(Opt(Style='ForeignEquity')),I1=1540, I2=1/90, 
                         sigma1=0.14,sigma2=0.18, g1=0.02,rho=-0.3,Type=c('Foreign', 'Domestic')){
  stopifnot(is.OptPx(o), o$Style$ForeignEquity, 
            is.numeric(I1), is.numeric(I2),is.numeric(g1),is.numeric(sigma1),
            is.numeric(sigma2), is.numeric(rho),is.character(Type))
  
  Type=match.arg(Type)
  isMax = switch(Type, 'Foreign'=TRUE, 'Domestic'=FALSE)
  
  d_f = (log(I1/o$K)+(o$r-g1-(sigma1^2)/2))/(sigma1*sqrt(o$ttm))
  d_1f = d_f+sigma1*(sqrt(o$ttm))
  A1 = I1*exp(-g1*o$ttm) * stats::pnorm(d_1f)
  A2 = -I1*exp(-g1*o$ttm) * stats::pnorm(-d_1f)
  B1 = o$K*exp(-o$r*o$ttm) * stats::pnorm(d_f)
  B2 = -o$K*exp(-o$r*o$ttm) * stats::pnorm(-d_f)
  C1 = I1*exp(rho*sigma1*sigma2-g1) * stats::pnorm(d_1f+rho*sigma2*sqrt(o$ttm))
  C2 = -I1*exp(rho*sigma1*sigma2-g1) * stats::pnorm(-d_1f-rho*sigma2*sqrt(o$ttm))
  D1 = o$K*exp(-o$r*o$ttm) * stats::pnorm(d_f+rho*sigma2*sqrt(o$ttm))
  D2 = -o$K*exp(-o$r*o$ttm) * stats::pnorm(-d_f-rho*sigma2*sqrt(o$ttm))
  # Calculate Price:
  
  if (Type=='Foreign' & o$Right$Name=='Call'){
    o$PxBS =A1-B1}
  
  if (Type =="Domestic"&  o$Right$Name=='Call'){
    o$PxBS =I2*(C1-D1)}
  
  if (Type == "Foreign"&  o$Right$Name=='Put'){
    o$PxBS =A2-B2}
  
  if (Type == "Domestic"& o$Right$Name=='Put'){
    o$PxBS =I2*(C2-D2)}
  return (o)
}


