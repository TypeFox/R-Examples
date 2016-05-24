#' Creates input for the function stochdom2 for stochastic dominance.
#' 
#' Stochastic dominance is a sophisticated comparison of two distributions of
#' stock market returns.  The dominating distribution is superior in terms of
#' mean, variance, skewness and kurtosis respectively, representing dominance
#' orders 1 to 4, without directly computing four moments.  Vinod(2008) sec. 4.3
#' explains the details.  The `wtdpapb' function creates the input
#' for stochdom2 which in turn computes the stochsastic dominance.
#' See Vinod (2004) for details about quantitative stochastic dominance.
#' 
#' @note  Internal function needed for stochastic dominance
#' 
#' @param xa {vector of (excess) returns for the first investment option A or
#' values of any random variable being compared to another.}
#' @param xb vector of returns for the second option B
#' @return 
#' \item{wpa}{weighted vector of probabilities for option A} 
#' \item{wpb}{weighted vector of probabilities for option B} 
#' \item{dj}{vector of distances when both values atre forced on the same support} 

#' @note In Vinod (2008) where the purpose of \code{wtdpapb} is to map from standard
#' `expected utility theory' weights to more sophisticated `non-expected utility
#' theory' weights using Prelec's (1998, Econometrica, p. 497) method.  These
#'  weights are not needed here. Hence we provide the function \code{prelec2}
#'  which does not use Prelec weights at all, thereby simplifying and speeding up
#'  the R code provided in Vinod (2008).
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Vinod, H. D.", "Hands-On Intermediate Econometrics 
#' Using R"  (2008) World Scientific Publishers: Hackensack, NJ.
#' \url{http://www.worldscibooks.com/economics/6895.html}
#'
#' @references Vinod, H. D. "Ranking Mutual Funds Using 
#' Unconventional Utility Theory and Stochastic Dominance,"
#' Journal of Empirical Finance Vol. 11(3) 2004, pp. 353-377.
#' 
#' @keywords stochastic dominance
#' @examples
#' 
#'  \dontrun{
#'  set.seed(234);x=sample(1:30);y=sample(5:34)
#'  wtdpapb(x,y)}
#' 
#' @export

wtdpapb <-
function(xa,xb){
#input: excess returns for mutual fund A & B
#output is a weighted pa and pb vectors (wpa, wpb) and dj vector
Ta=length(xa); Tb=length(xb);k=Ta+Tb
pa0=rep(1/Ta, Ta); pb0=rep(1/Tb, Tb)
xpapb=matrix(0,k,3)
for (i in 1:Ta){xpapb[i,1]=xa[i]
xpapb[i,2]=pa0[i]}
for (i in 1:Tb){xpapb[Ta+i,1]=xb[i]
xpapb[Ta+i,3]=pb0[i]}
pra=prelec2(n=Ta+Tb)
sm=sort_matrix(xpapb,1) #sort on first col of excess returns
pa=sm[,2]; pb=sm[,3]
wpa=pa*pra$pdif; wpb=pb*pra$pdif
dj=sm[,1]-sm[1,1]  #deviations from the minimum
#at sm[1,1] (sm is sorted)
list(wpa=wpa, wpb=wpb, dj=dj)}
