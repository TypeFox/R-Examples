###################################################################################################
#' Plot for Exact method of expected length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' The input can also be a range of values between 0 and 1.
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Confidence interval for \code{p} based on inverting equal-tailed
#' binomial tests with null hypothesis \eqn{H0: p = p0} using expected length of the \eqn{n + 1} intervals.
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;e=0.5;a=1;b=1
#' PlotexplEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=1;a=1;b=1 #Clopper-Pearson
#' PlotexplEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1);a=1;b=1 #Range including Mid-p and Clopper-Pearson
#' PlotexplEX(n,alp,e,a,b)
#' }
#' @export
##### 1.EXACT EMTHOD Expected Length for a given n and alpha level
PlotexplEX<-function(n,alp,e,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ewEX=NULL

  ELEX2=gexplEX(n,alp,e,a,b)
  ELEX2$e=as.factor(ELEX2$e)

  ggplot2::ggplot(ELEX2, ggplot2::aes(x=hp, y=ewEX, color=e))+
    ggplot2::labs(title = "Expected length of Exact method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 2)+
    ggplot2::labs(x = "p") +
    ggplot2::geom_line()

}

gexplEX<-function(n,alp,e,a,b)
{
nvar=length(e)

res <- data.frame()

for(i in 1:nvar)
{
  lu=ncf302(n,alp,e[i],a,b)
  res <- rbind(res,lu)
}
return(res)
}
ncf302<-function(n,alp,e,a,b)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LEX=0
UEX=0
s=5000
LEEX=0 								#LENGTH OF INTERVAL
ewiEX=matrix(0,k,s)						#Expected length quantity in sum
ewEX=0									#Expected Length

#EXACT METHOD
LEX[1]=0
UEX[1]= 1-((alp/(2*e))^(1/n))
LEX[k]=(alp/(2*e))^(1/n)
UEX[k]=1
LEEX[1]=1-((alp/(2*e))^(1/n))
LEEX[k]=1-((alp/(2*e))^(1/n))

for(i in 1:(k-2))
{
LEX[i+1]=exlim302l(x[i+1],n,alp,e)
UEX[i+1]=exlim302u(x[i+1],n,alp,e)
LEEX[i+1]=UEX[i+1]-LEX[i+1]
}
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiEX[i,j]=LEEX[i]*dbinom(i-1, n,hp[j])
}
ewEX[j]=sum(ewiEX[,j])						#Expected Length
}
#sumLEEX=sum(LEEX)
ELEX=data.frame(hp,ewEX,e)

return(ELEX)
}
#####TO FIND LOWER LIMITS
exlim302l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim302u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}

###############################################################################################################
#' Plot the Bayesian method of expected length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Plots of Bayesian Highest Probability Density (HPD) and two tailed
#' intervals using expected length of the \eqn{n + 1}
#' intervals for the Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1;a1=1;a2=1
#' PlotexplBA(n,alp,a,b,a1,a2)
#' }
#' @export
##### 8.BAYESIAN Expected Length for a given n and alpha level
PlotexplBA<-function(n,alp,a,b,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1)>1 || a1<0 ) stop("'a1' has to be greater than or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2)>1 || a2<0 ) stop("'a2' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LBAQ=0
UBAQ=0
LBAH=0
UBAH=0
s=5000
LEBAQ=0 								#LENGTH OF INTERVAL
LEBAH=0
ewiBAQ=matrix(0,k,s)						#Expected length quantity in sum
ewBAQ=0
ewiBAH=matrix(0,k,s)						#Expected length quantity in sum
ewBAH=0									#Expected Length

#library(TeachingDemos)				#To get HPDs
for(i in 1:k)
{
#Quantile Based Intervals
LBAQ[i]=qbeta(alp/2,x[i]+a1,n-x[i]+a2)
UBAQ[i]=qbeta(1-(alp/2),x[i]+a1,n-x[i]+a2)

LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[1]
UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[2]

LEBAQ[i]=UBAQ[i]-LBAQ[i]
LEBAH[i]=UBAH[i]-LBAH[i]
}
#sumLEBAQ=sum(LEBAQ)
#sumLEBAH=sum(LEBAH)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiBAQ[i,j]=LEBAQ[i]*dbinom(i-1, n,hp[j])
ewiBAH[i,j]=LEBAH[i]*dbinom(i-1, n,hp[j])

}
ewBAQ[j]=sum(ewiBAQ[,j])
ewBAH[j]=sum(ewiBAH[,j])						#Expected Length
}
ELBAQ=data.frame(hp,ew=ewBAQ,method="Quantile")
ELBAH=data.frame(hp,ew=ewBAH,method="HPD")

df.ba=rbind(ELBAQ,ELBAH)

ggplot2::ggplot(df.ba, ggplot2::aes(x=hp, y=ew))+
  ggplot2::labs(title = "Expected length of Bayesian Quantile & HPD based methods") +
  ggplot2::labs(y = "Expected Length") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(ggplot2::aes(color=method)) +
  ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 2)+
  ggplot2::scale_colour_manual(name='Heading',
                               values=c('Quantile' ='black',
                                        'HPD' = 'red'),
                               guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(1,1))))

}

#############################################################################################################
#' Plots the Expected length using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for the expected length of \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1;b=1;
#' PlotexplAll(n,alp,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplAll<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  #### Calling functions and creating df
  df.eall=  explAll(n,alp,a,b)

  ggplot2::ggplot(df.eall, ggplot2::aes(x=hp, y=ew))+
    ggplot2::labs(title = "Expected length of 6 base methods") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 1)

}

#############################################################################################################
#' Plots the expected length for Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type intervals using sum of length of the \eqn{n + 1}
#'  intervals
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' PlotexplWD(n,alp,a,b)
#' }
#' @export
##### 1.WALD sum of length for a given n and alpha level
PlotexplWD<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  df.wd=  gexplWD(n,alp,a,b)
  ddf.wd = lengthWD(n,alp,a,b)
  df.wd$gMean=ddf.wd$explMean
  df.wd$gMax=ddf.wd$explMax
  df.wd$gUL=ddf.wd$explMean+ddf.wd$explSD
  df.wd$gLL=ddf.wd$explMean-ddf.wd$explSD


  ggplot2::ggplot(data=df.wd, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Wald method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE)  +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))


}

#############################################################################################################
#' Plots the expected length for Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of score test approach using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' PlotexplSC(n,alp,a,b)
#' }
#' @export
##### 2.SCORE - sum of length for a given n and alpha level
PlotexplSC<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  df.sc=  gexplSC(n,alp,a,b)
  ddf.sc = lengthSC(n,alp,a,b)
  df.sc$gMean=ddf.sc$explMean
  df.sc$gMax=ddf.sc$explMax
  df.sc$gUL=ddf.sc$explMean+ddf.sc$explSD
  df.sc$gLL=ddf.sc$explMean-ddf.sc$explSD

  ggplot2::ggplot(data=df.sc, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Score method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

  }

#############################################################################################################
#' Plots ArcSine method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Wald-type interval for the arcsine transformation of the parameter
#' \code{p} using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' PlotexplAS(n,alp,a,b)
#' }
#' @export
##### 3. ARC SINE - sum of length for a given n and alpha level
PlotexplAS<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  df.as=  gexplAS(n,alp,a,b)

  ddf.as = lengthAS(n,alp,a,b)
  df.as$gMean=ddf.as$explMean
  df.as$gMax=ddf.as$explMax
  df.as$gUL=ddf.as$explMean+ddf.as$explSD
  df.as$gLL=ddf.as$explMean-ddf.as$explSD

  ggplot2::ggplot(data=df.as, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of ArcSine method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}

#############################################################################################################
#' Plots Logit Wald method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Wald-type interval based on the logit
#' transformation of \code{p} using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' PlotexplLT(n,alp,a,b)
#' }
#' @export
##### 4.LOGIT-WALD - sum of length for a given n and alpha level
PlotexplLT<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  df.lt=  gexplLT(n,alp,a,b)
  ddf.lt = lengthLT(n,alp,a,b)
  df.lt$gMean=ddf.lt$explMean
  df.lt$gMax=ddf.lt$explMax
  df.lt$gUL=ddf.lt$explMean+ddf.lt$explSD
  df.lt$gLL=ddf.lt$explMean-ddf.lt$explSD

  ggplot2::ggplot(data=df.lt, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Logit Wald method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}

#############################################################################################################
#' Plots Wald-T method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of approximate method based on a t_approximation of the
#' standardized point estimator using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' PlotexplTW(n,alp,a,b)
#' }
#' @export
##### 5.t-WALD - sum of length for a given n and alpha level
PlotexplTW<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  df.tw=  gexplTW(n,alp,a,b)
  ddf.tw = lengthTW(n,alp,a,b)
  df.tw$gMean=ddf.tw$explMean
  df.tw$gMax=ddf.tw$explMax
  df.tw$gUL=ddf.tw$explMean+ddf.tw$explSD
  df.tw$gLL=ddf.tw$explMean-ddf.tw$explSD

  ggplot2::ggplot(data=df.tw, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Wald-T method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))


}

#############################################################################################################
#' Plots likelihood Ratio method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Likelihood ratio limits using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' PlotexplLR(n,alp,a,b)
#' }
#' @export
#####6.LIKELIHOOD RATIO - sum of length for a given n and alpha level
PlotexplLR<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=NULL

  df.lr=  gexplLR(n,alp,a,b)
  ddf.lr = lengthLR(n,alp,a,b)
  df.lr$gMean=ddf.lr$explMean
  df.lr$gMax=ddf.lr$explMax
  df.lr$gUL=ddf.lr$explMean+ddf.lr$explSD
  df.lr$gLL=ddf.lr$explMean-ddf.lr$explSD

  ggplot2::ggplot(data=df.lr, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Likelihood Ratio method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}
