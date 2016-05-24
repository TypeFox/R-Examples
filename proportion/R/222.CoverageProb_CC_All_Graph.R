#############################################################################################################
#' Graphs  of Coverage Probability for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The graphs of the  Coverage Probability of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCAll(n,alp,c,a,b,t1,t2)
#' }
#' @export
##### 1.All methods - Coverage Probability 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
PlotcovpCAll<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL


  #### Calling functions and creating df
  df1    = gcovpCWD(n,alp,c,a,b,t1,t2)
  df2    = gcovpCAS(n,alp,c,a,b,t1,t2)
  df3    = gcovpCSC(n,alp,c,a,b,t1,t2)
  df4    = gcovpCLT(n,alp,c,a,b,t1,t2)
  df5    = gcovpCTW(n,alp,c,a,b,t1,t2)

  g.df= rbind(df1,df2,df3,df4,df5)
  g.df$t1=t1
  g.df$t2=t2
  g.df$alp=alp

  ggplot2::ggplot(g.df, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected methods") +
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
#' Graphs  of Coverage Probability for continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Wald method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCWD(n,alp,c,a,b,t1,t2)
#' }
#' @export
PlotcovpCWD<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  df1    = gcovpCWD(n,alp,c,a,b,t1,t2)
  WaldcovpA.df    = covpCWD(n,alp,c,a,b,t1,t2)
  df1$mcp=WaldcovpA.df$mcpCW
  df1$micp=WaldcovpA.df$micpCW
  df1$t1=t1
  df1$t2=t2
  df1$alp=alp

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Wald method") +
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
#' Graphs  of Coverage Probability for continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected ArcSine method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCAS(n,alp,c,a,b,t1,t2)
#' }
#' @export
PlotcovpCAS<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  df1    = gcovpCAS(n,alp,c,a,b,t1,t2)
  ArcSinecovpA.df = covpCAS(n,alp,c,a,b,t1,t2)

  df1$mcp=ArcSinecovpA.df$mcpCA
  df1$micp=ArcSinecovpA.df$micpCA
  df1$t1=t1
  df1$t2=t2
  df1$alp=alp

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected ArcSine method") +
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
#' Graphs  of Coverage Probability for continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Score method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCSC(n,alp,c,a,b,t1,t2)
#' }
#' @export
PlotcovpCSC<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  df1    = gcovpCSC(n,alp,c,a,b,t1,t2)
  ScorecovpA.df   = covpCSC(n,alp,c,a,b,t1,t2)

  df1$mcp=ScorecovpA.df$mcpCS
  df1$micp=ScorecovpA.df$micpCS
  df1$t1=t1
  df1$t2=t2
  df1$alp=alp

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Score method") +
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
#' Graphs  of Coverage Probability for continuity corrected Logistic Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Logistic Wald method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCLT(n,alp,c,a,b,t1,t2)
#' }
#' @export
PlotcovpCLT<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  df1    = gcovpCLT(n,alp,c,a,b,t1,t2)
  WaldLcovpA.df   = covpCLT(n,alp,c,a,b,t1,t2)

  df1$mcp=WaldLcovpA.df$mcpCLT
  df1$micp=WaldLcovpA.df$micpCLT
  df1$t1=t1
  df1$t2=t2
  df1$alp=alp

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Logistic Wald method") +
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
#' Graphs  of Coverage Probability for continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Wald-T method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCTW(n,alp,c,a,b,t1,t2)
#' }
#' @export
PlotcovpCTW<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  df1    = gcovpCTW(n,alp,c,a,b,t1,t2)

  AdWaldcovpA.df  = covpCTW(n,alp,c,a,b,t1,t2)

  df1$mcp=AdWaldcovpA.df$mcpCTW
  df1$micp=AdWaldcovpA.df$micpCTW
  df1$t1=t1
  df1$t2=t2
  df1$alp=alp

ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Wald-T method") +
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
