#' Calculates the sum of lengths for a specific LL and UL
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param hp - Hypothetical "p"
#' @details  Evaluation of intervals obtained from any method using sum
#' of the lengths for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length
#' @examples
#' n= 5;
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' hp=seq(0,1,by=0.01)
#' lengthGEN(n,LL,UL,hp)
#' @export
##### 1.Expected Length
lengthGEN<-function(n,LL,UL,hp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(hp)) stop("'hp' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if (any(hp>1) || any(hp<0)) stop("'hp' has to be between 0 and 1")

  ####INPUT n
  x=0:n
  k=n+1
  s=length(hp)
  ewi=matrix(0,k,s)						#Expected length quantity in sum
  ew=0									#Expected Length
  LE=0

  for(i in 1:k)
  {
    LE[i]=UL[i]-LL[i]
  }
  ####Expected Length

  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewi[i,j]=LE[i]*dbinom(i-1, n,hp[j])
    }
    ew[j]=sum(ewi[,j])						#Expected Length
  }

  sumLen=sum(LE)
  explMean=mean(ew)
  explSD=sd(ew)
  explMax=max(ew)
  explLL=explMean-(explSD)
  explUL=explMean+(explSD)
  df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
  return(df.Summary)
}

###########################################################################
#' Plots the sum of lengths for a specific LL and UL
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param hp - Hypothetical "p"
#' @details  Plot of intervals obtained from any method using sum
#' of the lengths for the \eqn{n + 1} intervals
#' @family Expected length
#' @examples
#' n= 5;
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' hp=seq(0,1,by=0.01)
#' PlotlengthGEN(n,LL,UL,hp)
#' @export
##### Plot Length
PlotlengthGEN<-function(n,LL,UL,hp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(hp)) stop("'hp' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if (any(hp>1) || any(hp<0)) stop("'hp' has to be between 0 and 1")
  method=sumLen=NULL

  full.df=lengthGEN(n,LL,UL,hp)
  full.df$method="General"

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values="orange") +
    ggplot2::labs(title = "Sum Length - General method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

###############################################################################################################
#' Sum of length  calculated using simulation method
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param s - Number of Hypothetical "p"
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The calculation of the sum of length for \code{n} given lower limit \code{LL} and  upper limit \code{UL}
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length
#' @examples
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' n= 5; s=5000; a=1; b=1;
#' lengthSIM(n,LL,UL,s,a,b)
#' @export
##### 2.Expected Length - Graph
lengthSIM<-function(n,LL,UL,s,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(s)) stop("'s' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if ((class(s) != "integer") & (class(s) != "numeric") || length(s)>1 || s<1  ) stop("'b' has to be greater than or equal to 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  ####INPUT n
  k=n+1
  LE=0
  ewi=matrix(0,k,s)						#sum of length quantity in sum
  ew=0

  for(i in 1:k)
  {
    LE[i]=UL[i]-LL[i]
  }
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewi[i,j]=LE[i]*dbinom(i-1, n,hp[j])
    }
    ew[j]=sum(ewi[,j])						#Expected Length
  }
  sumLen=sum(LE)
  explMean=mean(ew)
  explSD=sd(ew)
  explMax=max(ew)
  explLL=explMean-(explSD)
  explUL=explMean+(explSD)
  df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
  return(df.Summary)
}

###############################################################################################################
#' Plots the sum of length using calculated using simulation
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param s - Number of Hypothetical "p"
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The plot of the expected length for \code{n} given lower limit \code{LL} and  upper limit \code{UL}
#' @family Expected length
#' @examples
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' n= 5; s=5000; a=1; b=1;
#' PlotlengthSIM(n,LL,UL,s,a,b)
#' @export
##### 2.Expected Length - Graph
PlotlengthSIM<-function(n,LL,UL,s,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(s)) stop("'s' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if ((class(s) != "integer") & (class(s) != "numeric") || length(s)>1 || s<1  ) stop("'b' has to be greater than or equal to 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  sumLen=method=NULL

  ####INPUT n
  x=0:n
  k=n+1
  ewi=matrix(0,k,s)						#Expected length quantity in sum
  ew=0									#Expected Length
  LE=0

  for(i in 1:k)
  {
    LE[i]=UL[i]-LL[i]
  }
  sumLen=sum(LE)
  EL=data.frame(sumLen,method="Simulation")

  ggplot2::ggplot(EL, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values="orange") +
    ggplot2::labs(title = "Sum Length - Simulation method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

