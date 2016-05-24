#############################################################################################################
#' Plots the Expected length using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)  for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplAAll(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplAAll<-function(n,alp,h,a,b)
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
  df.new=  explAAll(n,alp,h,a,b)

  ggplot2::ggplot(df.new, ggplot2::aes(x=hp, y=ew))+
    ggplot2::labs(title = "Expected length of 6 adjusted methods") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 1)

}

#############################################################################################################
#' Plots the Expected length using adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplAWD(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplAWD<-function(n,alp,h,a,b)
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
  df.awd=  gexplAWD(n,alp,h,a,b)
  ddf.awd = lengthAWD(n,alp,h,a,b)
  df.awd$gMean=ddf.awd$explMean
  df.awd$gMax=ddf.awd$explMax
  df.awd$gUL=ddf.awd$explMean+ddf.awd$explSD
  df.awd$gLL=ddf.awd$explMean-ddf.awd$explSD

  ggplot2::ggplot(data=df.awd, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Wald method") +
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
#' Plots the Expected length using adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Score method
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplASC(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplASC<-function(n,alp,h,a,b)
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
  df.asc=  gexplASC(n,alp,h,a,b)
  ddf.asc = lengthASC(n,alp,h,a,b)
  df.asc$gMean=ddf.asc$explMean
  df.asc$gMax=ddf.asc$explMax
  df.asc$gUL=ddf.asc$explMean+ddf.asc$explSD
  df.asc$gLL=ddf.asc$explMean-ddf.asc$explSD

  ggplot2::ggplot(data=df.asc, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Score method") +
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
#' Plots the Expected length using adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplATW(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplATW<-function(n,alp,h,a,b)
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
  df.atw=  gexplATW(n,alp,h,a,b)
  ddf.atw = lengthATW(n,alp,h,a,b)
  df.atw$gMean=ddf.atw$explMean
  df.atw$gMax=ddf.atw$explMax
  df.atw$gUL=ddf.atw$explMean+ddf.atw$explSD
  df.atw$gLL=ddf.atw$explMean-ddf.atw$explSD

  ggplot2::ggplot(data=df.atw, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Wald-T method") +
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
#' Plots the Expected length using adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplALT(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplALT<-function(n,alp,h,a,b)
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
  df.alt=  gexplALT(n,alp,h,a,b)
  ddf.alt = lengthALT(n,alp,h,a,b)
  df.alt$gMean=ddf.alt$explMean
  df.alt$gMax=ddf.alt$explMax
  df.alt$gUL=ddf.alt$explMean+ddf.alt$explSD
  df.alt$gLL=ddf.alt$explMean-ddf.alt$explSD

  ggplot2::ggplot(data=df.alt, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Logit Wald method") +
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
#' Plots the Expected length using adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted ArcSine method
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplAAS(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplAAS<-function(n,alp,h,a,b)
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
  df.aas=  gexplAAS(n,alp,h,a,b)
  ddf.aas = lengthAAS(n,alp,h,a,b)
  df.aas$gMean=ddf.aas$explMean
  df.aas$gMax=ddf.aas$explMax
  df.aas$gUL=ddf.aas$explMean+ddf.aas$explSD
  df.aas$gLL=ddf.aas$explMean-ddf.aas$explSD

  ggplot2::ggplot(data=df.aas, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Logit Wald method") +
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
#' Plots the Expected length using adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Likelihood Ratio method
#' @family Expected length  of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplALR(n,alp,h,a,b)
#' }
#' @export
##### 9.All methods - Expected length
PlotexplALR<-function(n,alp,h,a,b)
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
  df.alr=  gexplALR(n,alp,h,a,b)
  ddf.alr = lengthALR(n,alp,h,a,b)
  df.alr$gMean=ddf.alr$explMean
  df.alr$gMax=ddf.alr$explMax
  df.alr$gUL=ddf.alr$explMean+ddf.alr$explSD
  df.alr$gLL=ddf.alr$explMean-ddf.alr$explSD

  ggplot2::ggplot(data=df.alr, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Likelihood Ratio method") +
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
