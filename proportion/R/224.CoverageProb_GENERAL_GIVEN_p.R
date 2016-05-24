#' General Coverage Probability given hypothetical "p" Coverage probability for CI
#' obtained from any method and discrete hypothetical p
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param alp - Alpha value (significance level required)
#' @param hp - Hypothetical "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of intervals obtained from any method using
#' coverage probability, root mean square statistic, and the proportion of proportion
#' lies within the desired level of coverage for the \eqn{n + 1} intervals and pre-defined
#' space for the parameter \code{p}
#' @return A dataframe with
#'  \item{mcp}{  Mean Coverage Probability}
#'  \item{micp }{   Minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Mean Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family General methods for coverage probability
#' @examples
#' \dontrun{
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' hp=seq(0,1,by=0.0001)
#' n= 5; alp=0.05; t1=0.93; t2=0.97
#' covpGEN(n,LL,UL,alp,hp,t1,t2)
#' }
#' @export
##### 1.Coverage Probability Related Outputs
covpGEN<-function(n,LL,UL,alp,hp,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(hp)) stop("'hp' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (any(hp>1) || any(hp<0) ) stop("'hp' has to be between 0 and 1")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
k=n+1
s=length(hp)
cp=matrix(0,k,s)
ct=matrix(0,k,s)							#Cover Pbty quantity in sum
cpp=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0
###CRITICAL VALUES
#cv=qnorm(1-(alp/2), mean = 0, sd = 1)
####COVERAGE PROBABILITIES
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LL[i] && hp[j] < UL[i])
{
cp[i,j]=dbinom(i-1, n,hp[j])
ct[i,j]=1
}
}
cpp[j]=sum(cp[,j])
RMSE_N1[j]=(cpp[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cpp[j]&&cpp[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CP=data.frame(hp,cpp)
mcp=mean(cpp)
micp=min(cpp)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cpp[j]-mcp)^2
RMSE_Mi1[j]=(cpp[j]-micp)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcp,micp,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Plot of simulation based Coverage Probability with discrete values for p
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param alp - Alpha value (significance level required)
#' @param hp - Hypothetical "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Graphical evaluation of intervals obtained from any method using
#' coverage probability, root mean square statistic, and the proportion of proportion
#' lies within the desired level of coverage for the \eqn{n + 1} intervals and pre-defined
#' space for the parameter \code{p}
#' @family General methods for coverage probability
#' @examples
#' \dontrun{
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' hp=seq(0,1,by=0.0001)
#' n= 5; alp=0.05; t1=0.93; t2=0.97
#' PlotcovpGEN(n,LL,UL,alp,hp,t1,t2)
#' }
#' @export
##### 2.  Coverage Probability - Graph
PlotcovpGEN<-function(n,LL,UL,alp,hp,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(hp)) stop("'hp' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (any(hp>1) || any(hp<0) ) stop("'hp' has to be between 0 and 1")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=cp=cpp=mcp=micp=NULL

  ####INPUT n
 # x=0:n
  k=n+1
  s=length(hp)
  cp=matrix(0,k,s)
  ct=matrix(0,k,s)							#Cover Pbty quantity in sum
  cpp=0								#Coverage probabilty
  ctr=0
  ###CRITICAL VALUES
  #cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  ####COVERAGE PROBABILITIES
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LL[i] && hp[j] < UL[i])
      {
        cp[i,j]=dbinom(i-1, n,hp[j])
        ct[i,j]=1
      }
    }
    cpp[j]=sum(cp[,j])
    if(t1<cpp[j]&&cpp[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CP=data.frame(hp,cpp)
  CP$mcp=mean(cpp)
  CP$micp=min(cpp)					#Mean Cov Prob


ggplot2::ggplot(CP, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(x = "p") +
  ggplot2::labs(title = "Coverage Probability - General method") +
  ggplot2::geom_line(ggplot2::aes(color="Coverage Probability"))+
  ggplot2::geom_point(ggplot2::aes(color="CP Values"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
  ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
  ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp),color="Confidence Level"),linetype = 2)+
  ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
  ggplot2::scale_colour_manual(name='Heading',
                        values=c('Coverage Probability'='red',
                                 'CP Values'='red',
                                 'Minimum Coverage'='black',
                                 'Mean Coverage'='blue',
                                 'Confidence Level'='brown'),
                        guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1,1,1,1),
                                                     shape=c(NA, NA, 16,NA,NA))))

}
