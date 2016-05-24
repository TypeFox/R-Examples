#############################################################################################################
#' Graphs of Coverage Probability - exact method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e -  Exact method indicator (1:Clop-Pear,0.5:MID-p)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  graphs of basic Coverage Probability methods
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; e=0.5; a=1;b=1; t1=0.93;t2=0.97 # Mid-p
#' PlotcovpEX(n,alp,e,a,b,t1,t2)
#' n= 10; alp=0.05; e=1; a=1;b=1; t1=0.93;t2=0.97 #Clop-Pear
#' PlotcovpEX(n,alp,e,a,b,t1,t2)
#' n=5; alp=0.05;
#' e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpEX(n,alp,e,a,b,t1,t2)
#' }
#' @export
##### 7. EXACT METHOD - Coverage Probability
PlotcovpEX=function(n,alp,e,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10 ) stop("Plot of only 10 intervals of 'e' is possible")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  method=Value=hp=cp=cpp=mcp=micp=NULL

  if(length(e)>1){
    dfex=gcovpEX(n,alp,e,a,b,t1,t2)
    exdf=dfex[,c(1,2,5)]
    exdf$e=as.factor(exdf$e)
    exdf$t1=t1
    exdf$t2=t2
    exdf$alp=alp

ggplot2::ggplot(exdf, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(title = "Coverage Probability for exact method for multiple e values") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
  ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
  ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
  ggplot2::geom_line(ggplot2::aes(color=e)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)
  }
  else{
    dfex=gcovpEX(n,alp,e,a,b,t1,t2)
dfex1=data.frame(micp=dfex$micpEX[1]	,mcp=dfex$mcpEX[1]	)
dfex1$alp=alp

ggplot2::ggplot(dfex, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(title = "Coverage Probability of exact method") +
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(ggplot2::aes(color="Coverage Probability"))+
  ggplot2::geom_point(ggplot2::aes(color="CP Values"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp),color="Confidence Level"),linetype = 2)+
  ggplot2::geom_hline(data=dfex1,ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
  ggplot2::geom_hline(data=dfex1,ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
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
}

#############################################################################################################
#' Graphs of  Coverage Probability of the Bayesian method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  The  graphs of  Coverage Probability of Bayesian method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97;a1=1;a2=1
#' PlotcovpBA(n,alp,a,b,t1,t2,a1,a2)
#' }
#' @export
#8.BAYESIAN
PlotcovpBA<-function(n,alp,a,b,t1,t2,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1)>1 || a1<0  ) stop("'a1' has to be greater than or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2)>1 || a2<0  ) stop("'a2' has to be greater than or equal to 0")
  ID=method=Value=hp=cp=cpp=mcpBAQ=micpBAQ=mcpBAH=micpBAH=NULL

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LBAQ=0
UBAQ=0
LBAH=0
UBAH=0
s=5000
cpBAQ=matrix(0,k,s)
ctBAQ=matrix(0,k,s)							#Cover Pbty quantity in sum
cppBAQ=0								#Coverage probabilty
ctr=0

cpBAH=matrix(0,k,s)
ctBAH=matrix(0,k,s)							#Cover Pbty quantity in sum
cppBAH=0								#Coverage probabilty
ctrH=0

##############
#library(TeachingDemos)				#To get HPDs
for(i in 1:k)
{
#Quantile Based Intervals
LBAQ[i]=qbeta(alp/2,x[i]+a1,n-x[i]+a2)
UBAQ[i]=qbeta(1-(alp/2),x[i]+a1,n-x[i]+a2)

LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[1]
UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[2]

}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LBAQ[i] && hp[j] < UBAQ[i])
{
cpBAQ[i,j]=dbinom(i-1, n,hp[j])
ctBAQ[i,j]=1
}
if(hp[j] > LBAH[i] && hp[j] < UBAH[i])
{
cpBAH[i,j]=dbinom(i-1, n,hp[j])
ctBAH[i,j]=1
}

}
cppBAQ[j]=sum(cpBAQ[,j])						#Coverage Probability
if(t1<cppBAQ[j]&&cppBAQ[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

cppBAH[j]=sum(cpBAH[,j])						#Coverage Probability
if(t1<cppBAH[j]&&cppBAH[j]<t2) ctrH=ctrH+1		#tolerance for cov prob - user defined

}
CPBAQ=data.frame(hp,cpp=cppBAQ,method="Quantile")
CPBAH=data.frame(hp,cpp=cppBAH,method="HPD")

df.new=rbind(CPBAQ,CPBAH)
df.new$mcpBAQ=mean(cppBAQ)
df.new$micpBAQ=min(cppBAQ)					#Mean Cov Prob

df.new$mcpBAH=mean(cppBAH)
df.new$micpBAH=min(cppBAH)					#Mean Cov Prob
df.new$alp=alp					#Mean Cov Prob


ggplot2::ggplot(df.new, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(title = "Coverage Probability of Bayesian methods") +
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(ggplot2::aes(color=method)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp),color="Confidence Level"),linetype = 2)+
  ggplot2::geom_hline(ggplot2::aes(yintercept=micpBAQ,color="Minimum Coverage Quantile"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=mcpBAQ,color="Mean Coverage Quantile"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=micpBAH,color="Minimum Coverage HPD"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=mcpBAH,color="Mean Coverage HPD"))+
  ggplot2::scale_colour_manual(name='Heading',
                               values=c('Quantile' ='black',
                                        'HPD' = 'red',
                                        'Minimum Coverage Quantile'='red',
                                        'Mean Coverage Quantile'='blue',
                                        'Minimum Coverage HPD'='black',
                                        'Mean Coverage HPD'='cyan',
                                        'Confidence Level'='brown'),
                               guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1,1,1,1,1,1))))


}

#############################################################################################################
#' Graphs of basic Coverage Probability 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  graphs of basic Coverage Probability methods
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpAll(n,alp,a,b,t1,t2)
#' }
#' @export
##### 9.  Coverage Probability - Graph
PlotcovpAll<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  df1=gcovpW(n,alp,a,b,t1,t2)
  df2=gcovpS(n,alp,a,b,t1,t2)
  df3=gcovpA(n,alp,a,b,t1,t2)
  df4=gcovpLT(n,alp,a,b,t1,t2)
  df5=gcovpTW(n,alp,a,b,t1,t2)
  df6=gcovpL(n,alp,a,b,t1,t2)

 nndf=  rbind(df1,df2,df3,df4,df5,df6)
 nndf$t1=t1
 nndf$t2=t2
 nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(title = "Coverage Probability for 6 base methods") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Wald method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpWD(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpWD<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  Waldcovp.df    = covpWD(n,alp,a,b,t1,t2)

  nndf=gcovpW(n,alp,a,b,t1,t2)
  nndf$mcp=Waldcovp.df$mcpW
  nndf$micp=Waldcovp.df$micpW
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base ArcSine method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpAS(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpAS<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ArcSinecovp.df = covpAS(n,alp,a,b,t1,t2)

  nndf=gcovpA(n,alp,a,b,t1,t2)
  nndf$mcp=ArcSinecovp.df$mcpA
  nndf$micp=ArcSinecovp.df$micpA
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for ArcSine method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Likelihood Ratio method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpLR(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpLR<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  LRcovp.df      = covpLR(n,alp,a,b,t1,t2)
 # ss1 = data.frame(method = LRcovp.df$method, MeanCP=LRcovp.df$mcpL, MinCP= LRcovp.df$micpL, RMSE_N=LRcovp.df$RMSE_N,RMSE_M=LRcovp.df$RMSE_M,RMSE_MI=LRcovp.df$RMSE_MI,tol=LRcovp.df$tol)

  nndf=gcovpL(n,alp,a,b,t1,t2)
  nndf$mcp=LRcovp.df$mcpL
  nndf$micp=LRcovp.df$micpL
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Likelihood Ratio method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Score method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpSC(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpSC<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  Scorecovp.df   = covpSC(n,alp,a,b,t1,t2)
  #ss1 = data.frame( MeanCP=Scorecovp.df$mcpS, MinCP= Scorecovp.df$micpS, RMSE_N=Scorecovp.df$RMSE_N,RMSE_M=Scorecovp.df$RMSE_M,RMSE_MI=Scorecovp.df$RMSE_MI,tol=Scorecovp.df$tol)

  nndf=gcovpS(n,alp,a,b,t1,t2)
  nndf$mcp=Scorecovp.df$mcpS
  nndf$micp=Scorecovp.df$micpS
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Score method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Logit Wald method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpLT(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpLT<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  WaldLcovp.df   = covpLT(n,alp,a,b,t1,t2)

  nndf=gcovpLT(n,alp,a,b,t1,t2)
  nndf$mcp=WaldLcovp.df$mcpLT
  nndf$micp=WaldLcovp.df$micpLT
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Logit Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Wald-T method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpTW(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpTW<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  AdWaldcovp.df  = covpTW(n,alp,a,b,t1,t2)

  nndf=gcovpTW(n,alp,a,b,t1,t2)
  nndf$mcp=AdWaldcovp.df$mcpTW
  nndf$micp=AdWaldcovp.df$micpTW
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Wald-T method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}
