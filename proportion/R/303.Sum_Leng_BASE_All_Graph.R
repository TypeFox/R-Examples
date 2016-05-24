#############################################################################################################
#' Plots the length summary for 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots for 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine).
#' Sum of the length is shown as bar graph.
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1;b=1;
#' PlotlengthAll(n,alp,a,b)
#' }
#' @export
PlotlengthAll<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthAll(n,alp,a,b)

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("gray", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length -6 Base methods") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the length summary for exact method
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
#' n=5; alp=0.05;e=0.5;a=1;b=1 # Mid-p
#' PlotlengthEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=1;a=1;b=1 #Clopper-Pearson
#' PlotlengthEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1);a=1;b=1 #Range including Mid-p and Clopper-Pearson
#' PlotlengthEX(n,alp,e,a,b)
#' }
#' @export
##### 1.EXACT EMTHOD Expected Length for a given n and alpha level
PlotlengthEX<-function(n,alp,e,a,b)
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
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthEX(n,alp,e,a,b)
  full.df$e=as.factor(full.df$e)

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = e, y = sumLen,  fill=e)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Exact method") +
    ggplot2::labs(x = "e") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

###############################################################################################################
#' Plot the Bayesian method of length summary
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Plots of Bayesian Highest Probability Density (HPD) and two tailed
#' intervals using expected length of the \eqn{n + 1}
#' intervals for the Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @return A dataframe with
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1;a1=1;a2=1
#' PlotlengthBA(n,alp,a,b,a1,a2)
#' @export
##### 8.BAYESIAN Expected Length for a given n and alpha level
PlotlengthBA<-function(n,alp,a,b,a1,a2)
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
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthBA(n,alp,a,b,a1,a2)

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c(  "#56B4E9", "orange")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Quantile & HPD of Bayesian method") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the length summary for Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type intervals using sum of length of the \eqn{n + 1}
#'  intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotlengthWD(n,alp,a,b)
#' @export
##### 1.WALD sum of length for a given n and alpha level
PlotlengthWD<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthWD(n,alp,a,b)
  full.df$method="Wald"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Wald method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the length summary for Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of score test approach using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotlengthSC(n,alp,a,b)
#' @export
##### 2.SCORE - sum of length for a given n and alpha level
PlotlengthSC<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthSC(n,alp,a,b)
  full.df$method="Score"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Score method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots ArcSine method of sum of length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Wald-type interval for the arcsine transformation of the parameter
#' \code{p} using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotlengthAS(n,alp,a,b)
#' @export
##### 3. ARC SINE - sum of length for a given n and alpha level
PlotlengthAS<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthAS(n,alp,a,b)
  full.df$method="ArcSine"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - ArcSine method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots Logit Wald method of sum of length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Wald-type interval based on the logit
#' transformation of \code{p} using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotlengthLT(n,alp,a,b)
#' @export
##### 4.LOGIT-WALD - sum of length for a given n and alpha level
PlotlengthLT<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthLT(n,alp,a,b)
  full.df$method="Logit Wald"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Logit Wald method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots Wald-T method of sum of length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of approximate method based on a t_approximation of the
#' standardized point estimator using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotlengthTW(n,alp,a,b)
#' @export
##### 5.t-WALD - sum of length for a given n and alpha level
PlotlengthTW<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthTW(n,alp,a,b)
  full.df$method="Wald-T"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Wald-T method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots likelihood Ratio method of sum of length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Likelihood ratio limits using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotlengthLR(n,alp,a,b)
#' @export
#####6.LIKELIHOOD RATIO - sum of length for a given n and alpha level
PlotlengthLR<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthLR(n,alp,a,b)
  full.df$method="Likelihood Ratio"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "green", "#56B4E9", "orange","#F0E442", "#CC79A7",
                 "cyan4", "pink", "cyan", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Likelihood Ratio method") +
    ggplot2::labs(x = "method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}
