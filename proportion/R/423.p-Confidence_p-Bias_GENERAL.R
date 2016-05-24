#' Performs p-Confidence and p-Bias estimation only using n,
#' lower limit and upper limit for general method
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @details  Evaluation of intervals obtained from any method using
#' p-confidence and p-bias for the  \eqn{n + 1} intervals
#' @return A dataframe with
#' \describe{
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' }
#' @family General methods for p-Confidence and p-Bias
#' @examples
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' n=5;
#' pCOpBIGEN(n,LL,UL)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
##### 1. p-confidence and p-bias
pCOpBIGEN<-function(n,LL,UL)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) < n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) < n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")

####INPUT n
#x=0:n
k=n+1
pcon=0					#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LL[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LL[i]))
pconC[i-1]=2*pbinom(i-1, n, UL[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=(1-max(pcon[i-1],pconC[i-1]))*100 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=max(0,pbia1[i-1])*100
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}
###############################################################################################
#' Plots p-Confidence and p-Bias estimation for general method
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @details  Plot of the general method - intervals obtained from any method using
#' p-confidence and p-bias for the  \eqn{n + 1} intervals
#' @family General methods for p-Confidence and p-Bias
#' @examples
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' n=5;
#' PlotpCOpBIGEN(n,LL,UL)
#' @export
##### 2. p-confidence and p-bias
PlotpCOpBIGEN<-function(n,LL,UL)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) < n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) < n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  x=Value=Heading=mark=NULL

  gdf=pCOpBIGEN(n,LL,UL)
  W1 = data.frame(x=gdf$x1, Value=gdf$pconf, Heading="pconf")
  W2 = data.frame(x=gdf$x1, Value=gdf$pbias, Heading="pbias")
  ndf=rbind(W1,W2)

  ggplot2::ggplot(ndf, ggplot2::aes(x=x, y=Value))+
    ggplot2::labs(title = "p-Confidence & p-Bias - General method") +
    ggplot2::facet_grid(Heading ~ .,scales="free_y") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=W1,ggplot2::aes(color="pConfidence "))+
    ggplot2::geom_point(data=W1, ggplot2::aes(color="pconf Values"))+
    ggplot2::geom_line(data=W2,ggplot2::aes(color="pbias"))+
    ggplot2::geom_point(data=W2,ggplot2::aes(color="pbias Values"))+
    ggplot2::scale_colour_manual(name='Heading',
                                 values=c('pConfidence '='red',
                                          'pconf Values'='red',
                                          'pbias'='blue',
                                          'pbias Values'='blue'),
                                 guide='Title') +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(1,1,1,1),
                                                                       shape=c(NA, 16,16,NA))))

}
###############################################################################################

