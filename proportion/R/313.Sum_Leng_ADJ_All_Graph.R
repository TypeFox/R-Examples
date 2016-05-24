#############################################################################################################
#' Plots the length summary for 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots for 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine).
#' Sum of the length is shown as bar graph.
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2; a=1;b=1;
#' PlotlengthAAll(n,alp,h,a,b)
#' }
#' @export
PlotlengthAAll<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthAAll(n,alp,h,a,b)

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("gray", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted methods") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the summary length using adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the summary length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotlengthAWD(n,alp,h,a,b)
#' @export
##### 9.All methods - summary length
PlotlengthAWD<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthAWD(n,alp,h,a,b)
  full.df$method="Adj-Wald"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted Wald") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the summary length using adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the summary length of adjusted Score method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotlengthASC(n,alp,h,a,b)
#' @export
##### 9.All methods - summary length
PlotlengthASC<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthAAS(n,alp,h,a,b)
  full.df$method="Adj-Score"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted Score") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the summary length using adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the summary length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotlengthATW(n,alp,h,a,b)
#' @export
##### 9.All methods - summary length
PlotlengthATW<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthATW(n,alp,h,a,b)
  full.df$method="Adj-Wald-T"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted Wald-T") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the summary length using adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the summary length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotlengthALT(n,alp,h,a,b)
#' @export
##### 9.All methods - summary length
PlotlengthALT<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthALT(n,alp,h,a,b)
  full.df$method="Adj-Logit Wald"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted Logit Wald") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the summary length using adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the summary length of adjusted ArcSine method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotlengthAAS(n,alp,h,a,b)
#' @export
##### 9.All methods - summary length
PlotlengthAAS<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthAAS(n,alp,h,a,b)
  full.df$method="Adj-ArcSine"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted ArcSine") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the summary length using adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the summary length of adjusted Likelihood Ratio method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotlengthALR(n,alp,h,a,b)
#' @export
##### 9.All methods - summary length
PlotlengthALR<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthALR(n,alp,h,a,b)
  full.df$method="Adj-Likelihood Ratio"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length - Adjusted Likelihood Ratio") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}
