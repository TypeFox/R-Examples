#' @title Rainbow option valuation via Black-Scholes (BS) model
#' @description Rainbow Option via Black-Scholes (BS) model
#' @details 
#' Two types of Rainbow options are priced: \code{'Max'} and \code{'Min'}. 
#' @author Chengwei Ge,Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}
#' @param Type Rainbow option type: 'Max' or 'Min'.   
#' @param S1 A spot price of the underlying security 1 (usually S1)
#' @param S2 A spot price of the underlying security 2 (usually S2)
#' @param sigma1 a vector of implied volatilities for the associated security 1
#' @param sigma2 a vector of implied volatilities for the associated security 2 
#' @param D1 A percent yield per annum from the underlying security 1
#' @param D2 A percent yield per annum from the underlying security 2
#' @param rho is the correlation between asset 1 and asset 2
#' @return A list of class \code{RainbowBS} consisting of the original \code{OptPx} object 
#' and the option pricing parameters \code{S1}, \code{Type}, \code{isMax}, and \code{isMin}
#' as well as the computed price \code{PxBS}.
#' 
#' @references Zhang Peter G., \emph{Exotic Options}, 2nd ed, 1998.  
#' @examples
#' (o = RainbowBS())$PxBS
#'   
#'   o = OptPx(Opt(Style = 'Rainbow',  Right = "Put"), r = 0.08)
#'   RainbowBS(o, S1=100, S2=95, D1=0,D2=0,sigma1=0.15,sigma2=0.2, rho=0.75,Type='Min')
#'   
#'   o = OptPx(Opt(Style = 'Rainbow', K = 102, ttm = 1, Right = "Put"), r = 0.08)
#'   RainbowBS(o, S1=100, S2=95, D1=0,D2=0,sigma1=0.15,sigma2=0.2, rho=0.75,Type='Min')
#'   
#'   o=OptPx(Opt(Style = 'Rainbow', K = 102, ttm = 1, Right = "Put"), r = 0.08)
#'   RainbowBS(o, S1=100, S2=95, D1=0,D2=0,sigma1=0.15,sigma2=0.2, rho=0.75,Type='Max')
#'   
#'   o=OptPx(Opt(Style = 'Rainbow', K = 102, ttm = 1, Right = "Call"), r = 0.08)
#'   RainbowBS(o, S1=100, S2=95, D1=0,D2=0,sigma1=0.15,sigma2=0.2, rho=0.75,Type='Min')
#'   
#'   o=OptPx(Opt(Style = 'Rainbow', K = 102, ttm = 1, Right = "Call"), r = 0.08)
#'   RainbowBS(o, S1=100, S2=95, D1=0,D2=0,sigma1=0.15,sigma2=0.2, rho=0.75,Type='Max') 
#' @export
#' 
RainbowBS=function(
  o=OptPx(Opt(Style='Rainbow')),
  S1=100,S2=95,D1=0,D2=0,sigma1=0.15,sigma2=0.2, rho=0.75,Type=c('Max', 'Min')){

  stopifnot(is.OptPx(o), o$Style$Rainbow, 
            is.numeric(S1),is.numeric(S2),is.numeric(D1),is.numeric(D2),is.numeric(sigma1),
            is.numeric(sigma2), is.numeric(rho),is.character(Type))
  
  Type=match.arg(Type)
  isMax = switch(Type, 'Max'=TRUE, 'Min'=FALSE)  
  
  sigmaA<-sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2)
  rho1<-(rho*sigma2-sigma1)/sigmaA
  rho2<-(rho*sigma1-sigma2)/sigmaA
  d1 <- (log(S1/o$K)+(o$r-D1-sigma1^2/2)*o$ttm)/(sigma1*sqrt(o$ttm))
  d11<- d1+sigma1*sqrt(o$ttm)
  d2 <- (log(S2/o$K)+(o$r-D2-sigma2^2/2)*o$ttm)/(sigma2*sqrt(o$ttm))
  d22<- d2+sigma2*sqrt(o$ttm)
  d12 <- (log(S2/S1)+(D1-D2-sigmaA^2/2)*o$ttm)/(sigmaA*sqrt(o$ttm))
  d21 <- (log(S1/S2)+(D1-D2-sigmaA^2/2)*o$ttm)/(sigmaA*sqrt(o$ttm))
  A1<-S1*exp(-D1*o$ttm)*(pbnorm(d11,d12,rho1))
  B1<-S2*exp(-D2*o$ttm)*(pbnorm(d22,d21,rho2))
  C1<-o$K*exp(-o$r*o$ttm)*(pbnorm(d1,d2,rho))
  A2<-S1*exp(-D1*o$ttm)*(pbnorm(d11,-d12,-rho1))
  B2<-S2*exp(-D2*o$ttm)*(pbnorm(d22,-d21,-rho2))
  C2<-o$K*exp(-o$r*o$ttm)*(1-pbnorm(-d1,-d2,rho))
  A3<-S1*exp(-D1*o$ttm)*(stats::pnorm(d12))
  B3<-S2*exp(-D2*o$ttm)*(stats::pnorm(d21))  
  A4<-S1*exp(-D1*o$ttm)*(stats::pnorm(-d12))
  B4<-S2*exp(-D2*o$ttm)*(stats::pnorm(-d21))
  # Calculate Price:
  
  if (o$Right$Name=='Call' & Type=='Min'){
    o$PxBS =A1+B1-C1}
  
  if (o$Right$Name=='Call'& Type== 'Max'){
    o$PxBS =A2+B2-C2}
  
  if (o$Right$Name=='Put'& Type=='Min'){
    o$PxBS=A1+B1-C1+o$K*exp(-o$r*o$ttm)-(A3+B3)}
  
  if (o$Right$Name=='Put'& Type=='Max'){
    o$PxBS=A2+B2-C2+o$K*exp(-o$r*o$ttm)-(A4+B4)}
  
  return(o)
}
