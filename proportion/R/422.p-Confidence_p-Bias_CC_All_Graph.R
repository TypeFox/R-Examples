###########################################################################################################3
#' Plots  p-confidence and p-bias for a given n and alpha level for 5 continuity corrected methods (Wald, Wald-T,  Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots of p-Confidence and p-Bias estimation of 5 continuity corrected methods (Wald, Wald-T,  Score, Logit-Wald, ArcSine)
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;c=1/(2*n)
#' PlotpCOpBICAll(n,alp,c)
#' }
#' @export
####5 CC methods plots of p-confidence and p-bias
PlotpCOpBICAll<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  x=Value=Heading=mark=val=NULL

  nAll = pCOpBICAll(n,alp,c)
  pc=data.frame(x=nAll$x1, val=nAll$pconf, Heading=nAll$method, mark="pconf",minmax=min(nAll$pbias, nAll$pconf))
  pb=data.frame(x=nAll$x1, val=nAll$pbias, Heading=nAll$method, mark="pbias",minmax=max(nAll$pbias, nAll$pconf))
  nndf=rbind(pc,pb)
  cdfWDc = subset(nndf, Heading == "CC-Wald" & mark== "pconf")
  cdfWDb = subset(nndf, Heading == "CC-Wald" & mark== "pbias")
  cdfSCc = subset(nndf, Heading == "CC-Score" & mark== "pconf")
  cdfSCb = subset(nndf, Heading == "CC-Score" & mark== "pbias")
  cdfASc = subset(nndf, Heading == "CC-ArcSine" & mark== "pconf")
  cdfASb = subset(nndf, Heading == "CC-ArcSine" & mark== "pbias")
  cdfLTc = subset(nndf, Heading == "CC-Logit-Wald" & mark== "pconf")
  cdfLTb = subset(nndf, Heading == "CC-Logit-Wald" & mark== "pbias")
  cdfTWc = subset(nndf, Heading == "CC-Wald-T" & mark== "pconf")
  cdfTWb = subset(nndf, Heading == "CC-Wald-T" & mark== "pbias")

  ggplot2::ggplot(nndf, ggplot2::aes(x=x, y=val))+
    ggplot2::labs(title = "p-Confidence & p-Bias - Continuity corrected methods") +
    ggplot2::facet_wrap(Heading ~ mark,scales="free_y",ncol=2) +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=cdfWDc,ggplot2::aes(color="pConf CC-Wald"))+
    ggplot2::geom_point(data=cdfWDc, ggplot2::aes(color="pConf CC-Wald Values"),shape=22)+
    ggplot2::geom_line(data=cdfWDb,ggplot2::aes(color="pbias CC-Wald"))+
    ggplot2::geom_point(data=cdfWDb,ggplot2::aes(color="pbias CC-Wald Values"),shape=23)+
    ggplot2::geom_line(data=cdfSCc,ggplot2::aes(color="pConf CC-Score"))+
    ggplot2::geom_point(data=cdfSCc, ggplot2::aes(color="pConf CC-Score Values"),shape=21)+
    ggplot2::geom_line(data=cdfSCb,ggplot2::aes(color="pbias CC-Score"))+
    ggplot2::geom_point(data=cdfSCb,ggplot2::aes(color="pbias CC-Score Values"),shape=16)+
    ggplot2::geom_line(data=cdfASc,ggplot2::aes(color="pConf CC-ArcSine"))+
    ggplot2::geom_point(data=cdfASc, ggplot2::aes(color="pConf CC-ArcSine Values"),shape=22)+
    ggplot2::geom_line(data=cdfASb,ggplot2::aes(color="pbias CC-ArcSine"))+
    ggplot2::geom_point(data=cdfASb,ggplot2::aes(color="pbias CC-ArcSine Values"),shape=23)+
    ggplot2::geom_line(data=cdfTWc,ggplot2::aes(color="pConf CC-Wald-T"))+
    ggplot2::geom_point(data=cdfTWc, ggplot2::aes(color="pConf CC-Wald-T Values"),shape=21)+
    ggplot2::geom_line(data=cdfTWb,ggplot2::aes(color="pbias CC-Wald-T"))+
    ggplot2::geom_point(data=cdfTWb,ggplot2::aes(color="pbias CC-Wald-T Values"),shape=16)+
    ggplot2::geom_line(data=cdfLTc,ggplot2::aes(color="pConf CC-Logit-Wald"))+
    ggplot2::geom_point(data=cdfLTc, ggplot2::aes(color="pConf CC-Logit-Wald Values"),shape=21)+
    ggplot2::geom_line(data=cdfLTb,ggplot2::aes(color="pbias CC-Logit-Wald"))+
    ggplot2::geom_point(data=cdfLTb,ggplot2::aes(color="pbias CC-Logit-Wald Values"),shape=16)+
    ggplot2::scale_colour_manual(name='Heading',
                                 values=c(
                                   'pConf CC-Wald'='red',
                                   'pConf CC-Wald Values'='red',
                                   'pbias CC-Wald'='black',
                                   'pbias CC-Wald Values'='black',
                                   'pConf CC-Score'='red',
                                   'pConf CC-Score Values'='red',
                                   'pbias CC-Score'='black',
                                   'pbias CC-Score Values'='black',
                                   'pConf CC-ArcSine'='red',
                                   'pConf CC-ArcSine Values'='red',
                                   'pbias CC-ArcSine'='black',
                                   'pbias CC-ArcSine Values'='black',
                                   'pConf CC-Logit-Wald'='red',
                                   'pConf CC-Logit-Wald Values'='red',
                                   'pbias CC-Logit-Wald'='black',
                                   'pbias CC-Logit-Wald Values'='black',
                                   'pConf CC-Wald-T'='red',
                                   'pConf CC-Wald-T Values'='red',
                                   'pbias CC-Wald-T'='black',
                                   'pbias CC-Wald-T Values'='black'))

}

######################################################################
#' Plots  p-confidence and p-bias for continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  p-confidence and p-bias plots for continuity corrected Wald method
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotpCOpBICWD(n,alp,c)
#' @export
PlotpCOpBICWD<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBICWD(n,alp,c)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - continuity corrected Wald method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  p-confidence and p-bias plots for continuity corrected ArcSine method
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotpCOpBICAS(n,alp,c)
#' @export
PlotpCOpBICAS<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBICAS(n,alp,c)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - continuity corrected ArcSine method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  p-confidence and p-bias plots for continuity corrected Wald-T method
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotpCOpBICTW(n,alp,c)
#' @export
PlotpCOpBICTW<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBICTW(n,alp,c)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - continuity corrected Wald-T method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  p-confidence and p-bias plots for continuity corrected Logit Wald method
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotpCOpBICLT(n,alp,c)
#' @export
PlotpCOpBICLT<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBICLT(n,alp,c)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - continuity corrected Logit Wald method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  p-confidence and p-bias plots for continuity corrected Score method
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' PlotpCOpBICSC(n,alp,c)
#' @export
PlotpCOpBICSC<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBICSC(n,alp,c)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - continuity corrected Score method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

