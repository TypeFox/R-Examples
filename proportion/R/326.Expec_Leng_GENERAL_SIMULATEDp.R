###############################################################################################################
#' Plots the expected length using calculated using simulation
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
#' PlotexplSIM(n,LL,UL,s,a,b)
#' @export
##### 2.Expected Length - Graph
PlotexplSIM<-function(n,LL,UL,s,a,b)
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
  hp=ew=method=explMean=explMax=explLL=explUL=NULL

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
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewi[i,j]=LE[i]*dbinom(i-1, n,hp[j])
}
ew[j]=sum(ewi[,j])						#Expected Length
}
explMean=mean(ew)
explSD=sd(ew)
explMax=max(ew)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
EL=data.frame(hp,ew,method="Simulation",explMean,explMax,explLL,explUL)

ggplot2::ggplot(data=EL, mapping=ggplot2::aes(x=hp, y=ew)) +
  ggplot2::labs(title = "Expected length of Simulation method") +
  ggplot2::labs(y = "Expected length") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show.legend  = TRUE)  +
  ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=explMean, fill="Mean"),color="orange"  ) +
  ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=explMax, fill="Max"),color="blue"  ) +
  ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=explLL, fill="Lower Limit"),color="cyan4"  ) +
  ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=explUL, fill="Upper Limit"),color="brown"  ) +
  ggplot2::scale_color_hue("Method") +
  ggplot2::scale_fill_manual(
    "Metric lines", values=c(1,1,1,1),
    guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
    labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}
