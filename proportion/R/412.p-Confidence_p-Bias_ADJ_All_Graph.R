#########################################################################################
#' Plots  p-confidence and p-bias for a given n and alpha level for
#' 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Plots of p-confidence and p-bias for 6 adjusted methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;h=2
#' PlotpCOpBIAAll(n,alp,h)
#' }
#' @export
####3.All methods plots of p-confidence and p-bias
PlotpCOpBIAAll<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  x=Value=Heading=val=mark=NULL

  nAll = pCOpBIAAll(n,alp,h)
  pc=data.frame(x=nAll$x1, val=nAll$pconf, Heading=nAll$method, mark="pconf",minmax=min(nAll$pbias, nAll$pconf))
  pb=data.frame(x=nAll$x1, val=nAll$pbias, Heading=nAll$method, mark="pbias",minmax=max(nAll$pbias, nAll$pconf))
  nndf=rbind(pc,pb)
  cdfWDc = subset(nndf, Heading == "Adj-Wald" & mark== "pconf")
  cdfWDb = subset(nndf, Heading == "Adj-Wald" & mark== "pbias")
  cdfSCc = subset(nndf, Heading == "Adj-Score" & mark== "pconf")
  cdfSCb = subset(nndf, Heading == "Adj-Score" & mark== "pbias")
  cdfASc = subset(nndf, Heading == "Adj-ArcSine" & mark== "pconf")
  cdfASb = subset(nndf, Heading == "Adj-ArcSine" & mark== "pbias")
  cdfLTc = subset(nndf, Heading == "Adj-Logit-Wald" & mark== "pconf")
  cdfLTb = subset(nndf, Heading == "Adj-Logit-Wald" & mark== "pbias")
  cdfLRc = subset(nndf, Heading == "Adj-Likelihood" & mark== "pconf")
  cdfLRb = subset(nndf, Heading == "Adj-Likelihood" & mark== "pbias")
  cdfTWc = subset(nndf, Heading == "Adj-Wald-T" & mark== "pconf")
  cdfTWb = subset(nndf, Heading == "Adj-Wald-T" & mark== "pbias")

  ggplot2::ggplot(nndf, ggplot2::aes(x=x, y=val))+
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted methods") +
    ggplot2::facet_wrap(Heading ~ mark,scales="free_y",ncol=2) +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=cdfWDc,ggplot2::aes(color="pConf Adj-Wald"))+
    ggplot2::geom_point(data=cdfWDc, ggplot2::aes(color="pConf Adj-Wald Values"),shape=22)+
    ggplot2::geom_line(data=cdfWDb,ggplot2::aes(color="pbias Adj-Wald"))+
    ggplot2::geom_point(data=cdfWDb,ggplot2::aes(color="pbias Adj-Wald Values"),shape=23)+
    ggplot2::geom_line(data=cdfSCc,ggplot2::aes(color="pConf Adj-Score"))+
    ggplot2::geom_point(data=cdfSCc, ggplot2::aes(color="pConf Adj-Score Values"),shape=21)+
    ggplot2::geom_line(data=cdfSCb,ggplot2::aes(color="pbias Adj-Score"))+
    ggplot2::geom_point(data=cdfSCb,ggplot2::aes(color="pbias Adj-Score Values"),shape=16)+
    ggplot2::geom_line(data=cdfASc,ggplot2::aes(color="pConf Adj-ArcSine"))+
    ggplot2::geom_point(data=cdfASc, ggplot2::aes(color="pConf Adj-ArcSine Values"),shape=22)+
    ggplot2::geom_line(data=cdfASb,ggplot2::aes(color="pbias Adj-ArcSine"))+
    ggplot2::geom_point(data=cdfASb,ggplot2::aes(color="pbias Adj-ArcSine Values"),shape=23)+
    ggplot2::geom_line(data=cdfTWc,ggplot2::aes(color="pConf Adj-Wald-T"))+
    ggplot2::geom_point(data=cdfTWc, ggplot2::aes(color="pConf Adj-Wald-T Values"),shape=21)+
    ggplot2::geom_line(data=cdfTWb,ggplot2::aes(color="pbias Adj-Wald-T"))+
    ggplot2::geom_point(data=cdfTWb,ggplot2::aes(color="pbias Adj-Wald-T Values"),shape=16)+
    ggplot2::geom_line(data=cdfLRc,ggplot2::aes(color="pConf Adj-Likelihood"))+
    ggplot2::geom_point(data=cdfLRc, ggplot2::aes(color="pConf Adj-Likelihood Values"),shape=22)+
    ggplot2::geom_line(data=cdfLRb,ggplot2::aes(color="pbias Adj-Likelihood"))+
    ggplot2::geom_point(data=cdfLRb,ggplot2::aes(color="pbias Adj-Likelihood Values"),shape=23)+
    ggplot2::geom_line(data=cdfLTc,ggplot2::aes(color="pConf Adj-Logit-Wald"))+
    ggplot2::geom_point(data=cdfLTc, ggplot2::aes(color="pConf Adj-Logit-Wald Values"),shape=21)+
    ggplot2::geom_line(data=cdfLTb,ggplot2::aes(color="pbias Adj-Logit-Wald"))+
    ggplot2::geom_point(data=cdfLTb,ggplot2::aes(color="pbias Adj-Logit-Wald Values"),shape=16)+
    ggplot2::scale_colour_manual(name='Heading',
                                 values=c(
                                   'pConf Adj-Wald'='red',
                                   'pConf Adj-Wald Values'='red',
                                   'pbias Adj-Wald'='black',
                                   'pbias Adj-Wald Values'='black',
                                   'pConf Adj-Score'='red',
                                   'pConf Adj-Score Values'='red',
                                   'pbias Adj-Score'='black',
                                   'pbias Adj-Score Values'='black',
                                   'pConf Adj-ArcSine'='red',
                                   'pConf Adj-ArcSine Values'='red',
                                   'pbias Adj-ArcSine'='black',
                                   'pbias Adj-ArcSine Values'='black',
                                   'pConf Adj-Logit-Wald'='red',
                                   'pConf Adj-Logit-Wald Values'='red',
                                   'pbias Adj-Logit-Wald'='black',
                                   'pbias Adj-Logit-Wald Values'='black',
                                   'pConf Adj-Wald-T'='red',
                                   'pConf Adj-Wald-T Values'='red',
                                   'pbias Adj-Wald-T'='black',
                                   'pbias Adj-Wald-T Values'='black',
                                   'pConf Adj-Likelihood'='red',
                                   'pConf Adj-Likelihood Values'='red',
                                   'pbias Adj-Likelihood'='black',
                                   'pbias Adj-Likelihood Values'='black'))

}

######################################################################
#' Plots  p-confidence and p-bias for adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  p-confidence and p-bias plots for adjusted Wald method
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' PlotpCOpBIAWD(n,alp,h)
#' @export
PlotpCOpBIAWD<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n)>1 || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIAWD(n,alp,h)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted Wald method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  p-confidence and p-bias plots for adjusted Likelihood Ratio method
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' PlotpCOpBIALR(n,alp,h)
#' @export
PlotpCOpBIALR<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n)>1 || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIALR(n,alp,h)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted Likelihood Ratio method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  p-confidence and p-bias plots for adjusted Wald-T method
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' PlotpCOpBIATW(n,alp,h)
#' @export
PlotpCOpBIATW<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n)>1 || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIATW(n,alp,h)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted Wald-T method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  p-confidence and p-bias plots for adjusted Logit Wald method
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' PlotpCOpBIALT(n,alp,h)
#' @export
PlotpCOpBIALT<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n)>1 || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIALT(n,alp,h)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted Logit Wald method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  p-confidence and p-bias plots for adjusted Score method
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' PlotpCOpBIASC(n,alp,h)
#' @export
PlotpCOpBIASC<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n)>1 || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIASC(n,alp,h)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted Score method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  p-confidence and p-bias plots for adjusted ArcSine method
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' PlotpCOpBIAAS(n,alp,h)
#' @export
PlotpCOpBIAAS<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n)>1 || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIAAS(n,alp,h)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Adjusted ArcSine method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}
