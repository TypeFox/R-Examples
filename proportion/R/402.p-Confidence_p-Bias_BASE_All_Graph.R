#################################################################################################
#' Plots of p-confidence and p-bias of Exact method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' The input can also be a range of values between 0 and 1.
#' @details  Evaluation of Confidence interval for \code{p} based on inverting equal-tailed
#' binomial tests with null hypothesis \eqn{H0: p = p0} using p-confidence and p-bias for
#' the \eqn{n + 1} intervals
#' @return A dataframe with
#' \describe{
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' }
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05;e=0.5; # Mid-p
#' PlotpCOpBIEX(n,alp,e)
#' n=5; alp=0.05;e=1; #Clopper-Pearson
#' PlotpCOpBIEX(n,alp,e)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1); #Range including Mid-p and Clopper-Pearson
#' PlotpCOpBIEX(n,alp,e)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 1.EXACT METHODS - p-confidence and p-bias for a given n and alpha level
PlotpCOpBIEX<-function(n,alp,e)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10 ) stop("Plot of only 10 intervals of 'e' is possible")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIEX(n,alp,e)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  b1=data.frame(e=CBEX$e)
  b2=data.frame(e=CBEX$e)
  bjoin=rbind(b1,b2)
  fg=rbind(W1,W2)
  gdf=cbind(fg,bjoin)
  gdf$e=as.factor(gdf$e)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value, colour = e)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Exact method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
      ggplot2::scale_colour_manual(values=c("blue", "red", "black", "cyan4", "deeppink",
                                            "orange","chartreuse4",
                                            "blueviolet" , "grey", "darksalmon", "tan1")) +
      ggplot2::geom_line(data=gdf,ggplot2::aes(color=e))

}

#######################################################################################################
#' Plots p-confidence and p-bias for Bayesian method given n and alpha level and priors a & b
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a1 - Shape parameter 1 for prior Beta distribution in Bayesian model
#' @param a2 - Shape parameter 2 for prior Beta distribution in Bayesian model
#' @details  Evaluation of Bayesian Highest Probability Density (HPD) and two tailed
#' intervals using p-confidence and p-bias for the \eqn{n + 1} intervals for the
#' Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05;a1=1;a2=1
#' PlotpCOpBIBA(n,alp,a1,a2)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
####2.BAYESIAN p-confidence and p-bias
PlotpCOpBIBA<-function(n,alp,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1)>1 || a1<0  ) stop("'a1' has to be greater than or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2)>1 || a2<0  ) stop("'a2' has to be greater than or equal to 0")
  x=Value=Heading=mark=NULL

  ndf=pCOpBIBA(n,alp,a1,a2)
  W1 = data.frame(x=ndf$x1, Value=ndf$pconfQ, Heading="pconfQ")
  W2 = data.frame(x=ndf$x1, Value=ndf$pbiasQ, Heading="pbiasQ")
  W3 = data.frame(x=ndf$x1, Value=ndf$pconfH, Heading="pconfH")
  W4 = data.frame(x=ndf$x1, Value=ndf$pbiasH, Heading="pbiasH")
  nBA=rbind(W1,W2,W3,W4)

  ggplot2::ggplot(nBA, ggplot2::aes(x=x, y=Value))+
    ggplot2::labs(title = "p-Confidence & p-Bias - Bayesian methods") +
    ggplot2::facet_grid(Heading ~ .,scales="free_y") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=W1,ggplot2::aes(color="pConf Quantile"))+
    ggplot2::geom_point(data=W1, ggplot2::aes(color="pconfQ Values"),shape=22)+
    ggplot2::geom_line(data=W2,ggplot2::aes(color="pbias Quantile"))+
    ggplot2::geom_point(data=W2,ggplot2::aes(color="pbiasQ Values"),shape=23)+
    ggplot2::geom_line(data=W3,ggplot2::aes(color="pConf HPD"))+
    ggplot2::geom_point(data=W3, ggplot2::aes(color="pconfHPD Values"),shape=21)+
    ggplot2::geom_line(data=W4,ggplot2::aes(color="pbias HPD"))+
    ggplot2::geom_point(data=W4,ggplot2::aes(color="pbiasHPD Values"),shape=16)+
    ggplot2::scale_colour_manual(name='Heading',
                                 values=c('pConf Quantile'='red',
                                          'pconfQ Values'='red',
                                          'pbias Quantile'='blue',
                                          'pbiasQ Values'='blue',
                                          'pConf HPD'='orange',
                                          'pconfHPD Values'='orange',
                                          'pbias HPD'='black',
                                          'pbiasHPD Values'='black'),
                                 guide='Title') +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(1,1,1,1,1,1,1,1),
                                                                       shape=c(NA, NA,16,23,NA,NA,21,22))))


}

###########################################################################################################3
#' Plots  p-confidence and p-bias for a given n and alpha level for
#' 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @family p-confidence and p-bias of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05
#' PlotpCOpBIAll(n,alp)
#' }
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
####3.All methods plots of p-confidence and p-bias
PlotpCOpBIAll<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=val=NULL

  nAll = pCOpBIAll(n,alp)
  pc=data.frame(x=nAll$x1, val=nAll$pconf, Heading=nAll$method, mark="pconf",minmax=min(nAll$pbias, nAll$pconf))
  pb=data.frame(x=nAll$x1, val=nAll$pbias, Heading=nAll$method, mark="pbias",minmax=max(nAll$pbias, nAll$pconf))
  nndf=rbind(pc,pb)
  cdfWDc = subset(nndf, Heading == "Wald" & mark== "pconf")
  cdfWDb = subset(nndf, Heading == "Wald" & mark== "pbias")
  cdfSCc = subset(nndf, Heading == "Score" & mark== "pconf")
  cdfSCb = subset(nndf, Heading == "Score" & mark== "pbias")
  cdfASc = subset(nndf, Heading == "ArcSine" & mark== "pconf")
  cdfASb = subset(nndf, Heading == "ArcSine" & mark== "pbias")
  cdfLTc = subset(nndf, Heading == "Logit-Wald" & mark== "pconf")
  cdfLTb = subset(nndf, Heading == "Logit-Wald" & mark== "pbias")
  cdfLRc = subset(nndf, Heading == "Likelihood" & mark== "pconf")
  cdfLRb = subset(nndf, Heading == "Likelihood" & mark== "pbias")
  cdfTWc = subset(nndf, Heading == "Wald-T" & mark== "pconf")
  cdfTWb = subset(nndf, Heading == "Wald-T" & mark== "pbias")

  ggplot2::ggplot(nndf, ggplot2::aes(x=x, y=val))+
    ggplot2::labs(title = "p-Confidence & p-Bias - All methods") +
    ggplot2::facet_wrap(Heading ~ mark,scales="free_y",ncol=2) +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=cdfWDc,ggplot2::aes(color="pConf Wald"))+
    ggplot2::geom_point(data=cdfWDc, ggplot2::aes(color="pConf Wald Values"),shape=22)+
    ggplot2::geom_line(data=cdfWDb,ggplot2::aes(color="pbias Wald"))+
    ggplot2::geom_point(data=cdfWDb,ggplot2::aes(color="pbias Wald Values"),shape=23)+
       ggplot2::geom_line(data=cdfSCc,ggplot2::aes(color="pConf Score"))+
    ggplot2::geom_point(data=cdfSCc, ggplot2::aes(color="pConf Score Values"),shape=21)+
    ggplot2::geom_line(data=cdfSCb,ggplot2::aes(color="pbias Score"))+
    ggplot2::geom_point(data=cdfSCb,ggplot2::aes(color="pbias Score Values"),shape=16)+
       ggplot2::geom_line(data=cdfASc,ggplot2::aes(color="pConf ArcSine"))+
    ggplot2::geom_point(data=cdfASc, ggplot2::aes(color="pConf ArcSine Values"),shape=22)+
    ggplot2::geom_line(data=cdfASb,ggplot2::aes(color="pbias ArcSine"))+
    ggplot2::geom_point(data=cdfASb,ggplot2::aes(color="pbias ArcSine Values"),shape=23)+
       ggplot2::geom_line(data=cdfTWc,ggplot2::aes(color="pConf Wald-T"))+
    ggplot2::geom_point(data=cdfTWc, ggplot2::aes(color="pConf Wald-T Values"),shape=21)+
    ggplot2::geom_line(data=cdfTWb,ggplot2::aes(color="pbias Wald-T"))+
    ggplot2::geom_point(data=cdfTWb,ggplot2::aes(color="pbias Wald-T Values"),shape=16)+
       ggplot2::geom_line(data=cdfLRc,ggplot2::aes(color="pConf Likelihood"))+
    ggplot2::geom_point(data=cdfLRc, ggplot2::aes(color="pConf Likelihood Values"),shape=22)+
    ggplot2::geom_line(data=cdfLRb,ggplot2::aes(color="pbias Likelihood"))+
    ggplot2::geom_point(data=cdfLRb,ggplot2::aes(color="pbias Likelihood Values"),shape=23)+
       ggplot2::geom_line(data=cdfLTc,ggplot2::aes(color="pConf Logit-Wald"))+
    ggplot2::geom_point(data=cdfLTc, ggplot2::aes(color="pConf Logit-Wald Values"),shape=21)+
    ggplot2::geom_line(data=cdfLTb,ggplot2::aes(color="pbias Logit-Wald"))+
    ggplot2::geom_point(data=cdfLTb,ggplot2::aes(color="pbias Logit-Wald Values"),shape=16)+
    ggplot2::scale_colour_manual(name='Heading',
                                 values=c(
                                    'pConf Wald'='red',
                                   'pConf Wald Values'='red',
                                   'pbias Wald'='black',
                                   'pbias Wald Values'='black',
                                    'pConf Score'='red',
                                   'pConf Score Values'='red',
                                   'pbias Score'='black',
                                   'pbias Score Values'='black',
                                    'pConf ArcSine'='red',
                                   'pConf ArcSine Values'='red',
                                   'pbias ArcSine'='black',
                                   'pbias ArcSine Values'='black',
                                    'pConf Logit-Wald'='red',
                                   'pConf Logit-Wald Values'='red',
                                   'pbias Logit-Wald'='black',
                                   'pbias Logit-Wald Values'='black',
                                    'pConf Wald-T'='red',
                                   'pConf Wald-T Values'='red',
                                   'pbias Wald-T'='black',
                                   'pbias Wald-T Values'='black',
                                'pConf Likelihood'='red',
                                          'pConf Likelihood Values'='red',
                                          'pbias Likelihood'='black',
                                          'pbias Likelihood Values'='black'))

}

######################################################################
#' Plots  p-confidence and p-bias for base Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for base Wald method
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' PlotpCOpBIWD(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
PlotpCOpBIWD<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIWD(n,alp)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Wald method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for base Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for base Score method
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' PlotpCOpBISC(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
PlotpCOpBISC<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBISC(n,alp)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Score method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for base ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for base ArcSine method
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' PlotpCOpBIAS(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
PlotpCOpBIAS<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBIAS(n,alp)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - ArcSine method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for base Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for base Logit Wald method
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' PlotpCOpBILT(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
PlotpCOpBILT<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBILT(n,alp)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Logit Wald method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for base Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for base Wald-T method
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' PlotpCOpBITW(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
PlotpCOpBITW<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBITW(n,alp)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Wald-T method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

######################################################################
#' Plots  p-confidence and p-bias for base Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  p-confidence and p-bias plots for base Likelihood Ratio method
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' PlotpCOpBILR(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
PlotpCOpBILR<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  x=Value=Heading=mark=NULL

  CBEX = pCOpBILR(n,alp)

  W1 = data.frame(x=CBEX$x1, Value=CBEX$pconf, Heading="pconf")
  W2 = data.frame(x=CBEX$x1, Value=CBEX$pbias, Heading="pbias")
  gdf=rbind(W1,W2)

  ggplot2::ggplot(gdf, ggplot2::aes(x=x, y=Value)) +
    ggplot2::facet_grid(Heading~.,scales="free_y") +
    ggplot2::labs(title = "p-Confidence & p-Bias - Likelihood Ratio method") +
    ggplot2::labs(y = "Confidence ") +
    ggplot2::labs(x = "No of successes") +
    ggplot2::geom_line(data=gdf,ggplot2::aes(x=x, y=Value))

}

