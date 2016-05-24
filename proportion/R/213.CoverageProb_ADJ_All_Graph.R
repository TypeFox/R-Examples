#############################################################################################################
#' Plots the Coverage Probability using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpAAll(n,alp,h,a,b,t1,t2)
#' }
#' @export
##### 9.All methods - Coverage Probability
PlotcovpAAll<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL


  #### Calling functions and creating df
  df.1    = gcovpAWD(n,alp,h,a,b,t1,t2)
  df.2    = gcovpAAS(n,alp,h,a,b,t1,t2)
  df.3    = gcovpALR(n,alp,h,a,b,t1,t2)
  df.4    = gcovpASC(n,alp,h,a,b,t1,t2)
  df.5    = gcovpALT(n,alp,h,a,b,t1,t2)
  df.6    = gcovpATW(n,alp,h,a,b,t1,t2)

  nndf=  rbind(df.1,df.2,df.3,df.4,df.5,df.6)
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of adjusted methods") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Wald method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpAWD(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpAWD<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  Waldcovp.df    = covpAWD(n,alp,h,a,b,t1,t2)
  nndf    = gcovpAWD(n,alp,h,a,b,t1,t2)
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  nndf$mcp=Waldcovp.df$mcpAW
  nndf$micp=Waldcovp.df$micpAW

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted ArcSine method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpAAS(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpAAS<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf    = gcovpAAS(n,alp,h,a,b,t1,t2)
  ArcSinecovp.df = covpAAS(n,alp,h,a,b,t1,t2)

  nndf$mcp=ArcSinecovp.df$mcpAA
  nndf$micp=ArcSinecovp.df$micpAA
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted ArcSine method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Likelihood Ratio method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpALR(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpALR<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf   = gcovpALR(n,alp,h,a,b,t1,t2)
  LRcovp.df      = covpALR(n,alp,h,a,b,t1,t2)

  nndf$mcp=LRcovp.df$mcpAL
  nndf$micp=LRcovp.df$micpAL
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Likelihood Ratio method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Score method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpASC(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpASC<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf   = gcovpASC(n,alp,h,a,b,t1,t2)
  Scorecovp.df   = covpASC(n,alp,h,a,b,t1,t2)

  nndf$mcp=Scorecovp.df$mcpAS
  nndf$micp=Scorecovp.df$micpAS
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Score method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Logistic Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Logistic Wald method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpALT(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpALT<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf    = gcovpALT(n,alp,h,a,b,t1,t2)
  WaldLcovp.df   = covpALT(n,alp,h,a,b,t1,t2)

  nndf$mcp=WaldLcovp.df$mcpALT
  nndf$micp=WaldLcovp.df$micpALT
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Logistic Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Wald-T method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpATW(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpATW<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf    = gcovpATW(n,alp,h,a,b,t1,t2)
  AdWaldcovp.df  = covpATW(n,alp,h,a,b,t1,t2)

  nndf$mcp=AdWaldcovp.df$mcpATW
  nndf$micp=AdWaldcovp.df$micpATW
  nndf$t1=t1
  nndf$t2=t2
  nndf$alp=alp

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Wald-T method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red") +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue") +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

