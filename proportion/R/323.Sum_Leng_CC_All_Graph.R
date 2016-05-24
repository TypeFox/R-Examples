#############################################################################################################
#' Plots the sum of length for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots for 5 continuity corrected  methods (Wald, Wald-T, Score, Logit-Wald, ArcSine).
#' Sum of the length is shown as bar graph.
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotlengthCAll(n,alp,c,a,b)
#' }
#' @export
PlotlengthCAll<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  full.df= lengthCAll(n,alp,c,a,b)

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("gray", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length  - Continuity corrected methods") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the sum of length using continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Wald method
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotlengthCWD(n,alp,c,a,b)
#' @export
##### 9.All methods - Expected length
PlotlengthCWD<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthCWD(n,alp,c,a,b)
  full.df$method="Continuity corrected Wald"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length  - Continuity corrected Wald") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the sum of length using continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Score method
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotlengthCSC(n,alp,c,a,b)
#' @export
##### 9.All methods - Expected length
PlotlengthCSC<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthCSC(n,alp,c,a,b)
  full.df$method="Continuity corrected Score"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length  - Continuity corrected Score") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the sum of length using continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected ArcSine method
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotlengthCAS(n,alp,c,a,b)
#' @export
##### 9.All methods - Expected length
PlotlengthCAS<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthCAS(n,alp,c,a,b)
  full.df$method="Continuity corrected ArcSine"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length  - Continuity corrected ArcSine") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the sum of length using continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Logit Wald method
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotlengthCLT(n,alp,c,a,b)
#' @export
##### 9.All methods - Expected length
PlotlengthCLT<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthCLT(n,alp,c,a,b)
  full.df$method="Continuity corrected Logit Wald"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length  - Continuity corrected Logit Wald") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

#############################################################################################################
#' Plots the sum of length using continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Wald-T method
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotlengthCTW(n,alp,c,a,b)
#' @export
##### 9.All methods - Expected length
PlotlengthCTW<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  hp=ew=method=gMean=gMax=gLL=gUL=explUL=explLL=sumLen=NULL

  #### Calling functions and creating df
  full.df= lengthCTW(n,alp,c,a,b)
  full.df$method="Continuity corrected Wald-T"

  g <- ggplot2::guide_legend("Mean")
  limits <- ggplot2::aes(ymax =explUL, ymin=explLL)
  cbPalette <- c("orange", "red", "#56B4E9", "orange","#F0E442", "#CC79A7")

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = sumLen,  fill=method)) +
    ggplot2::geom_bar(stat="identity",width=0.5) +
    ggplot2::scale_fill_manual(values=cbPalette) +
    ggplot2::labs(title = "Sum Length  - Continuity corrected Wald-T") +
    ggplot2::labs(x = "Method") +
    ggplot2::labs(y = "Sum of Length") +
    ggplot2::guides(colour = g, size = g, shape = g) +
    ggplot2::theme(legend.title = ggplot2::element_text(colour="black", size=11, face="bold"))+
    ggplot2::theme_classic()

}

