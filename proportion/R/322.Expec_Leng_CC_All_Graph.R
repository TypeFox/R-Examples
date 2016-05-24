#############################################################################################################
#' Plots the Expected length using  5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotexplCAll(n,alp,c,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplCAll<-function(n,alp,c,a,b)
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
  df.new    = explCAll(n,alp,c,a,b)

  ggplot2::ggplot(df.new, ggplot2::aes(x=hp, y=ew))+
    ggplot2::labs(title = "Expected length of continuity corrected methods") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 1)

}

#############################################################################################################
#' Plots the Expected length using continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Wald method
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotexplCWD(n,alp,c,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplCWD<-function(n,alp,c,a,b)
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
  df.cwd    = gexplCWD(n,alp,c,a,b)
  ddf.cwd = lengthCWD(n,alp,c,a,b)
  df.cwd$gMean=ddf.cwd$explMean
  df.cwd$gMax=ddf.cwd$explMax
  df.cwd$gUL=ddf.cwd$explMean+ddf.cwd$explSD
  df.cwd$gLL=ddf.cwd$explMean-ddf.cwd$explSD

  ggplot2::ggplot(data=df.cwd, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of continuity corrected Wald") +
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
#' Plots the Expected length using continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Score method
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotexplCSC(n,alp,c,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplCSC<-function(n,alp,c,a,b)
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
  df.csc    = gexplCSC(n,alp,c,a,b)
  ddf.csc = lengthCSC(n,alp,c,a,b)
  df.csc$gMean=ddf.csc$explMean
  df.csc$gMax=ddf.csc$explMax
  df.csc$gUL=ddf.csc$explMean+ddf.csc$explSD
  df.csc$gLL=ddf.csc$explMean-ddf.csc$explSD

  ggplot2::ggplot(data=df.csc, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of continuity corrected Score") +
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
#' Plots the Expected length using continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected ArcSine method
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotexplCAS(n,alp,c,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplCAS<-function(n,alp,c,a,b)
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
  df.cas    = gexplCAS(n,alp,c,a,b)
  ddf.cas = lengthCAS(n,alp,c,a,b)
  df.cas$gMean=ddf.cas$explMean
  df.cas$gMax=ddf.cas$explMax
  df.cas$gUL=ddf.cas$explMean+ddf.cas$explSD
  df.cas$gLL=ddf.cas$explMean-ddf.cas$explSD

  ggplot2::ggplot(data=df.cas, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of continuity corrected ArcSine") +
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
#' Plots the Expected length using continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Logit Wald method
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotexplCLT(n,alp,c,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplCLT<-function(n,alp,c,a,b)
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
  df.clt    = gexplCLT(n,alp,c,a,b)
  ddf.clt = lengthCLT(n,alp,c,a,b)
  df.clt$gMean=ddf.clt$explMean
  df.clt$gMax=ddf.clt$explMax
  df.clt$gUL=ddf.clt$explMean+ddf.clt$explSD
  df.clt$gLL=ddf.clt$explMean-ddf.clt$explSD

  ggplot2::ggplot(data=df.clt, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of continuity corrected Logit Wald") +
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
#' Plots the Expected length using continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correctection
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the expected length of continuity corrected Wald-T method
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' PlotexplCTW(n,alp,c,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplCTW<-function(n,alp,c,a,b)
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
  df.ctw    = gexplCTW(n,alp,c,a,b)
  ddf.ctw = lengthCTW(n,alp,c,a,b)
  df.ctw$gMean=ddf.ctw$explMean
  df.ctw$gMax=ddf.ctw$explMax
  df.ctw$gUL=ddf.ctw$explMean+ddf.ctw$explSD
  df.ctw$gLL=ddf.ctw$explMean-ddf.ctw$explSD

  ggplot2::ggplot(data=df.ctw, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of continuity corrected Wald-T") +
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
