#' @name autoplotTen
#' @title Generate a \code{ggplot} for a \code{survfit} or \code{ten} object
#'  
#' @include ten.R
#' @include print.R
#' @include autoplotTAP.R
#' @include sf.R
#' @include ci.R
#' @include nc.R
#'
#' @param object An object of class \code{survfit}, \code{ten} or \code{stratTen}.
#' @param ... Additional arguments (not implemented).
#' @param title Title for survival plot.
#' @param type \code{type="single"} (the default) plots single lines.
#'  \describe{
#'  \item{\code{type="CI"}}{Adds lines indicating
#'        confidence intervals (taken from \code{upper} and \code{lower}
#'        values of \code{survfit} object).
#'        \cr
#'        Higher values of \code{alpha} (opacity) are recommended for this,
#'        e.g. \code{alpha=0.8}.}
#'  \item{\code{type="fill"}}{Adds filled rectangles from the survival lines to
#'        the confidence intervals above.}
#' }
#' @param alpha Opacity of lines indicating confidence intervals
#' or filled rectangles. Should be in range \eqn{0-1}. Lower = more transparent.
#'  \cr
#' Larger values e.g. \code{alpha=0.7} are recommended for confidence
#' intervals.
#' @param ciLine Confidence interval line type. See 'line type specification' in 
#'  \cr
#' ?graphics::par
#' @param censShape Shape of marks to indicate censored onservations.
#'  \cr 
#' Default is \code{3} which gives vertical ticks.
#'  \cr 
#' Use \code{censShape=10} for circular marks. See 
#'  \cr
#' ?graphics::points
#' @param palette Options are taken from
#' \href{http://colorbrewer2.org/}{color_brewer}.
#'   \itemize{
#'     \item \code{palette="Dark2"} (the default) is recommended for
#'           \code{single} or \code{CI} plots.
#'     \item \code{palette="Set2"} is recommended for \code{type="fill"} plots.
#' }
#' @param jitter By default, \code{jitter="none"}.
#' \itemize{
#'  \item If \code{jitter="noEvents"}, adds some random, positive noise
#'   to survival lines with no events (i.e. all observations censored).
#'   This will bring them just above 1 on the y-axis, making them easier to see separately.
#'  \item If \code{jitter="all"} add some vertical 
#'   and horizontal noise to all survival lines. This can prevent overlapping 
#'   of lines for censoring.
#' }
#' @param tabTitle Table title.
#'  \cr \cr
#' \bold{--Axis arguments:}
#' @param xLab Label for \eqn{x} axis on survival plot.
#' @param timeTicks  Numbers to mark on the \eqn{x} axis of 
#' the survival plot and the table.
#' \describe{
#'  \item{\code{"major"}}{ 
#'   (the default) only the major \eqn{x}-axis (time) marks from the
#'   survival plot are are labelled on the plot and table.}
#'  \item{\code{"minor"}}{minor axis marks are labelled instead.}
#'  \item{\code{"days"}}{scale is \eqn{0, 7, 14, ..., t_{max}}}
#'  \item{\code{"months"}}{scale is \eqn{0, 12,, 24, ..., t_{max}}}
#'  \item{\code{"custom"}}{scale is given by \code{times} below}
#' }
#' @param times Vector of custom times to use for \eqn{x} axis.
#' @param yLab Label for \eqn{y} axis on survival plot.
#' @param yScale Display for point on \eqn{y} axis:
#' \describe{
#'  \item{\code{"perc"}}{Displays as percentages.}
#'  \item{\code{"frac"}}{Displays as fractions e.g. \eqn{0, 0.1, 0.2, ..., 1.0.}}
#' }
#' \bold{--Legend arguments:}
#'  \cr
#' @param legend If \code{legend=FALSE}, no legends will be produced
#' for the plot or table.
#' @param legTitle Legend title.
#' @param legLabs Legend labels. These can be used to replace the names
#' of the covariate groups ('strata' in the case of a \code{survfit} object).
#'  \cr
#' Should be given in the same order as those strata.
#' @param legOrd Legend order.
#'  \cr \cr
#' \bold{--Size arguments:}
#'  \cr
#'  Size arguments are passed to \code{ggplot2::element_text(size=)}.
#' @param titleSize Title size for survival plot.
#' @param  axisTitleSize Title size for axes.
#' @param  axisLabSize Title size for labels on axes.
#' @param survLineSize Survival line size.
#' 
#' @param legTitleSize Title size for legend.
#' @param legLabSize Legend labels width and height.
#' @param censSize Size of marks to indicate censored onservations.
#' @param fillLineSize Line size surrouding filled boxes.
#' @param tabTitleSize Table title text size.
#' @param tabLabSize Table legend text size.
#' @param nRiskSize Number at risk - text size.
#'  \cr \cr
#' \bold{--Arguments for autoplot.survfit only:}
#'  \cr
#' @param pVal If \code{pVal=TRUE}, adds \eqn{p} value from
#' log-rank test to plot
#' @param sigP No. of significant digits to display in \eqn{p} value.
#' Typically \eqn{1} to \eqn{3}.
#' @param pX Location of \eqn{p} value on \eqn{x} axis.
#'  \cr
#' Should be in the range of \eqn{0 - 1},
#' where value is to be placed relative to the maximum observed
#' time.
#'  \cr
#' E.g. \code{pX = 0.5} will place it half-way along \eqn{x}-axis
#' @param pY Location of \eqn{p} value on \eqn{y} axis.
#'  \cr
#' Should be in the range of \eqn{0 - 1}, as above.
#' 
#' @author Chris Dardis. \code{autoplot.survfit} based on existing work by
#' R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
#' 
#' @note \code{autoplot.survfit} may be deprecated after packageVersion 0.6.
#' Please try to use \code{autoplot.ten} instead.

#' 
#' @keywords hplot
#' @keywords survival
#' 
#' @seealso ?ggplot2::ggplot_build
#'
#' @rdname autoplotTen
#' @export
#' 
autoplot <- function (object, ...) UseMethod("autoplot")
#' 
#' @rdname autoplotTen
#' @method autoplot ten
#' @aliases autoplot.ten
#' @export
#' @examples
#' ## examples are slow to run; see vignette for output from these
#' \dontrun{
#' ### autoplot.ten
#' data("kidney", package="KMsurv")
#' t1 <- ten(survfit(Surv(time, delta) ~ type, data=kidney))
#' autoplot(t1)
#' autoplot(t1, type="fill", survLineSize=2, jitter="all")
#' autoplot(t1, timeTicks="months", 
#'  type="CI", jitter="all",
#'  legLabs=c("surgical", "percutaneous"),
#'  title="Time to infection following catheter placement \n
#'    by type of catheter, for dialysis patients",
#'  titleSize=10, censSize=2)$plot
#' t2 <- ten(survfit(Surv(time=time, event=delta) ~ 1, data=kidney))
#' autoplot(t2, legLabs="")$plot
#' autoplot(t2, legend=FALSE)
#' data("rectum.dat", package="km.ci")
#' t3 <- ten(survfit(Surv(time, status) ~ 1, data=rectum.dat))
#' ## change confidence intervals to log Equal-Precision confidence bands
#' ci(t3, how="nair", tL=1, tU=40)
#' autoplot(t3, type="fill", legend=FALSE)$plot
#' ## manually changing the output
#' t4 <- ten(survfit(Surv(time, delta) ~ type, data=kidney))
#' (a4 <- autoplot(t4, type="CI", alpha=0.8, survLineSize=2)$plot)
#' ## change default colors
#' a4 + list(ggplot2::scale_color_manual(values=c("red", "blue")),
#'           ggplot2::scale_fill_manual(values=c("red", "blue")))
#' ## change limits of y-axis
#' suppressMessages(a4 + ggplot2::scale_y_continuous(limits=c(0, 1)))
#' }
autoplot.ten <- function(object,
                         ...,
                         title="Marks show times with censoring",
                         type=c("single", "CI", "fill"),
                         alpha=0.05,
                         ciLine=10,
                         censShape=3,
                         palette=c("Dark2", "Set2", "Accent", "Paired",
                                   "Pastel1", "Pastel2", "Set1", "Set3"),
                         jitter=c("none", "noEvents", "all"),
                         tabTitle="Number at risk by time",
                         xLab="Time",
                         timeTicks=c("major", "minor", "days", "months", "custom"),
                         times=NULL,
                         yLab="Survival",
                         yScale=c("perc", "frac"),
                         legend=TRUE,
                         legTitle="Group",
                         legLabs=NULL,
                         legOrd=NULL,
                         titleSize=15,
                         axisTitleSize=15,
                         axisLabSize=10,
                         survLineSize=0.5,
                         censSize=5,
                         legTitleSize=10,
                         legLabSize=10,
                         fillLineSize=0.05,
                         tabTitleSize=15,
                         tabLabSize=5,
                         nRiskSize=5) {
    stopifnot(inherits(object, "ten"))
    stopifnot(alpha > 0 & alpha < 1)
    ## add no. censored
    nc(object)
    ## confidence intervals for plot
    dt1 <- data.table::copy(ci(object))
    dt1[, c("Sv", "SCV") := NULL]
    dt1 <- merge(object[, list(cg, t, nc)],
                 dt1,
                 all.x=FALSE,
                 all.y=FALSE,
                 by=c("cg", "t"))
    if (!is.null(legOrd)) {
        stopifnot(length(unique(legOrd))==length(unique(dt1[, cg])))
        stopifnot(all(legOrd %in% dt1[, seq.int(length(cg))]))
    }
    ## make two extra rows for each covariate group
    ## for t=0 to t=time of first event
    dt2 <- data.table::rbindlist(list(dt1[, .SD[1, ], by=cg],
                                      dt1[, .SD[1, ], by=cg]))
    ## set surv, upper and lower to one
    dt2[, c("S", "lower", "upper") := list(1), by=cg]
    ## set initial time and no. censored to zero
    dt2[seq.int(unique(dt2$cg)), c("t", "nc") := list(0L)]
    ## reorder to allow binding
    dt1 <- data.table::rbindlist(list(dt2, dt1))
### jitter
    jitter <- match.arg(jitter)
    ## for groups with no events add random no.to survival (by strata)
    if (jitter=="noEvents") {
        ## add column to indicate no. events by group
        dt1[, s1 := sum(n), by=list(cg)]
        dt1[s1==0, S := S + (stats::runif(1, 0.01, 0.05)), by=cg]
    }
    if(jitter=="all"){
        ## for groups with no events add random no.to survival (by strata)
        dt1[, S := S + (stats::runif(1, 0.01, 0.05)), by=cg]
        dt1[, t := abs(jitter(t, factor=0.5))]
    }
### 
    if (attr(object, "abbNames")) {
        na1 <- attr(object, "longNames")[, id]
        ## abbreviate function (used later)
        abbFn <- identity
    } else {
        na1 <- attr(object, "longNames")[, longName]
        abbFn <- as.integer
    }
    if (is.null(legLabs)) {
        dt1[, "cg" := factor(cg, labels=na1)]
    } else {
        stopifnot(length(legLabs)==length(unique(object$cg)))
        dt1[, "cg" := factor(cg, labels=legLabs)]
    }
    if (is.null(legOrd)) legOrd <- dt1[, seq.int(levels(cg))]
### 
### plot single lines only
### 
    g1 <- ggplot(data=dt1, aes(group=cg, color=cg, fill=cg)) +
        geom_step(aes(x=t, y=S), direction="hv", size=survLineSize)
###
    type <- match.arg(type)
    if (type=="CI") {
        g1 <- g1 +
            geom_step(aes(x=t, y=upper),
                      direction="hv", linetype=ciLine, alpha=alpha) +
            geom_step(aes(x=t, y=lower),
                      direction="hv", linetype=ciLine, alpha=alpha)
    }
    if (type=="fill") {
        ## copy dt1 to work allow further work
        dt2 <- dt1[, list(l=unique(lower),
                          u=unique(upper),
                          minT=as.numeric(min(t)),
                          t=as.numeric(t)
                          ), by=list(S, cg)]
        ## make max. time column
        dt2[, "maxT" := c(minT[2:length(minT)], Inf), by=cg]
        ## merge columns
        dt1 <- merge(dt1, dt2, by=c("t", "S", "cg"), all.y=TRUE)
        dt1 <- dt1[order(cg)]
        ## add shading
        g1 <- g1 + geom_rect(data=dt1, aes(x=NULL, y=NULL,
                                           ymax=S, ymin=l,
                                           xmax=maxT, xmin=minT,
                                           color=cg, group=cg, fill=cg),
                             alpha=alpha, size=fillLineSize) +
            geom_rect(data=dt1, aes(x=NULL, y=NULL,
                                    ymax=u, ymin=S,
                                    xmax=maxT, xmin=minT,
                                    color=cg, group=cg, fill=cg),
                      alpha=alpha, size=fillLineSize)
    }
    ## add lines to show times where subjects censored
    if (any(dt1[, nc >= 1])) {
        g1 <- g1 + geom_point(data=dt1[nc >= 1, ],
                              aes(x=t, y=S),
                              shape=censShape, size=censSize)
    }
### palette + legend
    ## use palette Dark2 for prominent shades
    ## (suitable for colorblind)
    ## use palette Set2 for lighter shades as large fill area
    palette <- match.arg(palette)
    g1 <- g1 + scale_color_brewer(type="qual",
                                  breaks=dt1[, levels(cg)[legOrd]],
                                  palette=palette,
                                  guide=guide_legend(
                                      title=legTitle)) +
        scale_fill_brewer(type="qual",
                          breaks=dt1[, levels(cg)[legOrd]],
                          palette=palette,
                          guide=guide_legend(
                              title=legTitle))
### scales
    g1 <- g1 + ggtitle(title)
    yScale <- match.arg(yScale)
    if (yScale=="frac") {
        g1 <- g1 + scale_y_continuous(yLab)
    } else {
        y1 <- ggplot_build(g1)$panel$ranges[[1L]]$y.major_source
        g1 <- g1 + scale_y_continuous(yLab,
                                      breaks=y1,
                                      labels=paste0(y1 * 100, "%"))
    }
    ## times to show
    ## use marks from existing plot
    timeTicks <- match.arg(timeTicks)
    x1 <- get("range", envir=get("range", envir=layer_scales(g1)$x))
    times1 <- switch(EXPR=timeTicks,
                     major=ggplot_build(g1)$panel$ranges[[1]]$x.major_source,
                     minor=ggplot_build(g1)$panel$ranges[[1]]$x.minor_source,
                     custom=NaN,
                     days=seq(from=min(x1), to=max(x1), by=7L),
                     months=seq(from=min(x1), to=max(x1), by=12L))
    if (is.nan(times1[1])) times1 <- times
    ## x axis
    g1 <- g1 +
        scale_x_continuous(name=xLab,
                           breaks=times1)
    ## font sizes
    g1 <- g1 +
        theme(title=element_text(size=titleSize),
              legend.text=element_text(size=legLabSize),
              legend.title=element_text(size=legTitleSize),
              axis.text=element_text(size=axisLabSize),
              axis.title=element_text(size=axisTitleSize))
### data for table of number at risk
    dt3 <- data.table::data.table("t"=times1)
    cg1 <- seq.int(attr(object, "ncg"))
    ## time, no. at risk, covariate group
    tnc1 <- lapply(cg1, FUN=function(cg1) {
        r1 <- data.table::setkey(object[abbFn(cg)==cg1, ncg, by=t], t)
        r1[dt3, roll=-Inf][, ncg]
    })
    tnc1 <- data.table::data.table(
        "t"=rep(times1, attr(object, "ncg")),
        "n"=unlist(tnc1),
        "cg"=as.factor(rep(na1, each=length(times1))))
    ## table
    g2 <- ggplot(data=tnc1, aes(x=t, y=cg, shape=cg)) +
        geom_point(size=0) +
        geom_text(aes(label=n), color=1, size=nRiskSize) +
        scale_x_continuous(name=xLab,
                           limits=c(0, max(dt1[, t])),
                           breaks=times1) +
        scale_y_discrete(name=legTitle, 
                         breaks=levels(tnc1$cg),
                         labels=levels(tnc1$cg)) +
        ggtitle(tabTitle) +
        theme(axis.text=element_text(size=axisLabSize),
              axis.title=element_text(size=axisTitleSize),
              plot.title=element_text(size=tabTitleSize),
              legend.title=element_text(size=tabLabSize),
              legend.text=element_text(size=tabLabSize)) +
        guides(shape=guide_legend(title=legTitle,
                                  keywidht=tabLabSize,
                                  keyheight=tabLabSize))
    ## remove legend
    if (!legend) {
        g1 <- g1 + theme(legend.position="none")
        g2 <- g2 + theme(legend.position="none")
    }
    res1 <- list("table"=g2,
                 "plot"=g1)
    class(res1) <- c("tableAndPlot", "list")
    return(res1)
}
#' 
#' @rdname autoplotTen
#' @method autoplot stratTen
#' @aliases autoplot.stratTen
#' @export
#' @examples
#' \dontrun{
#' data("pbc", package="survival")
#' t1 <- ten(Surv(time, status==2) ~ trt + strata(edema), data=pbc, abbNames=FALSE)
#' autoplot(t1)
#' }
autoplot.stratTen <- function(object,
                              ...,
                              title=NULL,
                              type=c("single", "CI", "fill"),
                              alpha=0.05,
                              ciLine=10,
                              censShape=3,
                              palette=c("Dark2", "Set2", "Accent", "Paired",
                                        "Pastel1", "Pastel2", "Set1", "Set3"),
                              jitter=c("none", "noEvents", "all"),
                              tabTitle="Number at risk by time",
                              xLab="Time",
                              timeTicks=c("major", "minor", "days", "months", "custom"),
                              times=NULL,
                              yLab="Survival",
                              yScale=c("perc", "frac"),
                              legend=TRUE,
                              legTitle="Group",
                              legLabs=NULL,
                              legOrd=NULL,
                              titleSize=15,
                              axisTitleSize=15,
                              axisLabSize=10,
                              survLineSize=0.5,
                              censSize=5,
                              legTitleSize=10,
                              legLabSize=10,
                              fillLineSize=0.05,
                              tabTitleSize=15,
                              tabLabSize=5,
                              nRiskSize=5) {
    res1 <- lapply(object,
                   autoplot,
                   title=title,
                   type=type,
                   alpha=alpha,
                   ciLine=ciLine,
                   censShape=3,
                   palette=palette,
                   jitter=jitter,
                   tabTitle=tabTitle,
                   xLab=xLab,
                   timeTicks=timeTicks,
                   times=times,
                   yLab=yLab,
                   yScale=yScale,
                   legend=TRUE,
                   legTitle=legTitle,
                   legLabs=legLabs,
                   legOrd=legOrd,
                   titleSize=titleSize,
                   axisTitleSize=axisTitleSize,
                   axisLabSize=axisLabSize,
                   survLineSize=survLineSize,
                   censSize=censSize,
                   legTitleSize=legTitleSize,
                   legLabSize=legLabSize,
                   fillLineSize=fillLineSize,
                   tabTitleSize=tabTitleSize,
                   tabLabSize=tabLabSize,
                   nRiskSize=nRiskSize)
    if (is.null(title)) {
        if (attr(object, "abbNames")) {
            title <- attr(object, "longNames")[, id]
        } else {
            title <- attr(object, "longNames")[, longName]
        }
    } else {
        title <- rep(title, length(object))
    }
    for (i in seq.int(length(object))){
        res1[[i]][[2]] <- res1[[i]][[2]] + ggplot2::ggtitle(title[i])
    }
    data.table::setattr(res1, "class", c("stratTableAndPlot", class(res1)))
    return(res1)
}
#'
#' @rdname autoplotTen
#' @method autoplot survfit
#' @aliases autoplot.survfit
#' @export
#' @examples
#' ### autoplot.survfit
#' \dontrun{
#' data(kidney, package="KMsurv")
#' s1 <- survfit(Surv(time, delta) ~ type, data=kidney)
#' autoplot(s1, type="fill", survLineSize=2)
#' autoplot(s1, type="CI", pVal=TRUE, pX=0.3,
#'  legLabs=c("surgical", "percutaneous"),
#'  title="Time to infection following catheter placement \n
#'    by type of catheter, for dialysis patients")$plot
#' s1 <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
#' autoplot(s1, legLabs="")$plot
#' autoplot(s1, legend=FALSE)$plot
#' data(rectum.dat, package="km.ci")
#' s1 <- survfit(Surv(time, status) ~ 1, data=rectum.dat)
#' ## change confidence intervals to log Equal-Precision confidence bands
#' if (require("km.ci")) {
#'  km.ci::km.ci(s1, method="logep")
#'  autoplot(s1, type="fill", legend=FALSE)$plot
#' }
#' ## manually changing the output
#' s1 <- survfit(Surv(time, delta) ~ type, data=kidney)
#' g1 <- autoplot(s1, type="CI", alpha=0.8, survLineSize=2)$plot
#' ## change default colors
#' g1 + ggplot2::scale_colour_manual(values=c("red", "blue")) +
#'     ggplot2::scale_fill_manual(values=c("red", "blue"))
#' ## change limits of y-axis
#' g1 + ggplot2::scale_y_continuous(limits=c(0, 1))
#' }
autoplot.survfit <- function(object,
                             ...,
                             title="Marks show times with censoring",
                             type=c("single", "CI", "fill"),
                             alpha=0.05,
                             ciLine=10,
                             censShape=3,
                             palette=c("Dark2", "Set2", "Accent", "Paired",
                                       "Pastel1", "Pastel2", "Set1", "Set3"),
                             jitter=c("none", "noEvents", "all"),
                             tabTitle="Number at risk by time",
                             xLab="Time",
                             timeTicks=c("major", "minor", "weeks", "months", "custom"),
                             times=NULL,
                             yLab="Survival",
                             yScale=c("perc", "frac"),
                             legend=TRUE,
                             legLabs=NULL,
                             legOrd=NULL,
                             legTitle="Group",
                             titleSize=15,
                             axisTitleSize=15,
                             axisLabSize=10,
                             survLineSize=0.5,
                             censSize=5,
                             legTitleSize=10,
                             legLabSize=10,
                             fillLineSize=0.05,
                             tabTitleSize=15,
                             tabLabSize=5,
                             nRiskSize=5,
                             pVal=FALSE,
                             sigP=1,
                             pX=0.1,
                             pY=0.1) {
    stopifnot(inherits(object, "survfit"))
    if (!is.null(legLabs) &! length(object$strata)==0){
        stopifnot(length(legLabs)==length(object$strata))
    }
    ## change names for strata to legLabs if required
    if (is.null(legLabs)) {
        stNames <- names(object$strata)
    } else {
        stNames <- legLabs
    }
    ## if only one strata (intercept only model)
    if (is.null(object$strata)) {
        if (is.null(legLabs)) {
            st1 <- as.factor(rep(1, length(object$time)))
        } else {
            stopifnot(length(legLabs)==1)
            st1 <- as.factor(rep(legLabs, length(object$time)))
        }
    } else {
        ## add vector for one strata according to number of rows of strata
        st1 <- unlist(sapply(1:length(object$strata),
                             function (i) rep(stNames[i], object$strata[i])))
    }
    ## create data.table with data from survfit
    ## add column for strata
    ## (using data.table here as avoids duplication when adding rows later)
    ## also rename strata as 'st' to avoid calling survival::function
    dt1 <- data.table::data.table(time=object$time,
                                  n.risk=object$n.risk,
                                  n.event=object$n.event,
                                  n.censor=object$n.censor,
                                  surv=object$surv,
                                  upper=object$upper,
                                  lower=object$lower,
                                  st=as.factor(st1))
    ## make two rows for each stratum
    ## for time=0 to time=time of first event
    dt2 <- data.table::rbindlist(list(dt1[, .SD[1, ], by=st],
                                      dt1[, .SD[1, ], by=st]))
    ## set n.event and n.censored to zero
    dt2[, c("n.event", "n.censor") := list(0), by=st]
    ## set surv, upper and lower to one
    dt2[, c("surv", "upper", "lower") := list(1), by=st]
    ## set first time to zero
    dt2[seq(length(unique(dt2$st))), "time" := (0L) ]
    ## reorder to allow binding
    data.table::setcolorder(dt2, names(dt1))
    dt1 <- data.table::rbindlist(list(dt2, dt1))
    if (is.null(legOrd)) legOrd <- dt1[, seq.int(levels(st))]
    ## 
    ## jitter
    ## 
    jitter <- match.arg(jitter)
    ## for groups with no events add random no.to survival (by strata)
    if (jitter=="noEvents") {
        ## add column to indicate no. events by group
        dt1[, s1 := sum(n.event), by=list(st)]
        dt1[s1==0, surv := surv+(runif(1, 0.01, 0.05)), by=st]
    }
    if(jitter=="all"){
        ## for groups with no events add random no.to survival (by strata)
        dt1[, surv := surv+(runif(1, 0.01, 0.05)), by=st]
    }
    ##
    dt1 <- dt1[order(st)]
    ## 
    ## plot single lines only
    ## 
    g1 <- ggplot(data=dt1, aes(group=st, colour=st, fill=st)) +
        geom_step(aes(x=time, y=surv), direction="hv", size=survLineSize)
    ##
    type <- match.arg(type)
    if (type=="CI"){
        g1 <- g1 +
            geom_step(aes(x=time, y=upper),
                      direction="hv", linetype=ciLine, alpha=alpha) +
            geom_step(aes(x=time, y=lower),
                      direction="hv", linetype=ciLine, alpha=alpha)
    }
    if (type=="fill"){
        ## copy dt1 to work allow further work
        dt2 <- dt1[, list(l=unique(lower),
                          u=unique(upper),
                          minT=as.numeric(min(time)),
                          time=as.numeric(time)
                          ), by=list(surv, st)]
        ## make max. time column
        dt2[, "maxT" := c(minT[2:length(minT)], NA), by=st]
        ## merge columns
        dt1 <- merge(dt1, dt2, by=c("time", "surv", "st"), all.y=TRUE)
        dt1 <- dt1[order(st)]
        ## add shading
        g1 <- g1 + geom_rect(data=dt1, aes(x=NULL, y=NULL,
                                           ymax=surv, ymin=l,
                                           xmax=maxT, xmin=minT,
                                           colour=st, group=st, fill=st),
                             alpha=alpha, size=fillLineSize) +
            geom_rect(data=dt1, aes(x=NULL, y=NULL,
                                    ymax=u, ymin=surv,
                                    xmax=maxT, xmin=minT,
                                    colour=st, group=st, fill=st),
                      alpha=alpha, size=fillLineSize)
    }
    ## add lines to show times where subjects censored
    if (any(dt1$n.censor >= 1)) {
        g1 <- g1 + geom_point(data=dt1[n.censor>=1, ],
                              aes(x=time, y=surv),
                              shape=censShape, size=censSize)
    }
    ## palette
    ## use palette Dark2 for prominent shades
    ## (suitable for colorblind)
    ## use palette Set2 for lighter shades as large fill area
    palette <- match.arg(palette)
    g1 <- g1 + scale_color_brewer(type="qual",
                                  breaks=dt1[, levels(cg)[legOrd]],
                                  palette=palette,
                                  guide=guide_legend(
                                      title=legTitle))
    g1 <- g1 + scale_fill_brewer(type="qual",
                                 breaks=dt1[, levels(cg)[legOrd]],
                                 palette=palette,
                                 guide=guide_legend(
                                     title=legTitle))
    ## scales
    g1 <- g1 + ggtitle(title)
    yScale <- match.arg(yScale)
    if (yScale=="frac") {
        g1 <- g1 + scale_y_continuous(yLab)
    } else {
        y1 <- ggplot_build(g1)$panel$ranges[[1L]]$y.major_source
        g1 <- g1 + scale_y_continuous(yLab,
                                      breaks=y1,
                                      labels=paste0(y1 * 100, "%"))
    }
    ## times to show
    timeTicks <- match.arg(timeTicks)
    x1 <- get("range", envir=get("range", envir=layer_scales(g1)$x))
    times1 <- switch(EXPR=timeTicks,
                     major=ggplot_build(g1)$panel$ranges[[1]]$x.major_source,
                     minor=ggplot_build(g1)$panel$ranges[[1]]$x.minor_source,
                     custom=NaN,
                     weeks=seq(from=min(x1), to=max(x1), by=7L),
                     months=seq(from=min(x1), to=max(x1), by=12L))
    if (is.nan(times1[1])) times1 <- times
    ## x axis
    g1 <- g1 +
        scale_x_continuous(name=xLab,
                           breaks=times1)
    ## font sizes
    g1 <- g1 +
        theme(title=element_text(size=titleSize),
              legend.text=element_text(size=legLabSize),
              legend.title=element_text(size=legTitleSize),
              axis.text = element_text(size = axisLabSize),
              axis.title = element_text(size = axisTitleSize))
    ## remove legend if required
    if(!legend) g1 <- g1 + theme(legend.position="none")
    ## p value for log-rank test (only if >=2 groups)
    if (pVal & !is.null(object$strata)) {
        sd1 <- survival::survdiff(eval(object$call$formula),
                                  data=eval(object$call$data))
        p1 <- stats::pchisq(sd1$chisq,
                            length(sd1$n) - 1,
                            lower.tail=FALSE)
        p1txt <- ifelse(p1 < 0.0001,
                        "Log-rank test \n p < 0.0001",
                        paste("Log-rank test \n p =", signif(p1, sigP)))
        g1 <- g1 + annotate(geom="text", 
                            x=pX * max(dt1$time),
                            y=pY,
                            label=p1txt,
                            size=legLabSize)
    }
    ## data for table
    dt3 <- data.table::data.table(
        time=summary(object, times = times1, extend = TRUE)$time,
        n.risk=summary(object, times = times1, extend = TRUE)$n.risk)
    ## if intercept-only model
    if (is.null(object$strata)) {
        dt3[, "st" := as.factor(rep(1, length(times1)))]
    } else {
        dt3[, "st" := summary(object, times=times1, extend=TRUE)$strata]
    }
    ## change names of strata to legend labels
    if(!is.null(legLabs)) dt3[, "st" := factor(st, labels=legLabs) ]
    ## table
    ## reverse here to plot in same order as in main plot
    g2 <- ggplot(data=dt3, aes(x=time, y=st, shape=st)) +
        geom_point(size=0) +
        geom_text(aes(label=n.risk), colour=1, size=nRiskSize) +
        scale_x_continuous(name=xLab,
                           limits=c(0, max(object$time)),
                           breaks=times1) +
        ## reverse here to plot in same order as in main plot
        scale_y_discrete(name=legTitle,
                         breaks=levels(dt3$st),
                         labels=levels(dt3$st)) +
        ggtitle(tabTitle) +
        theme(axis.text = element_text(size=axisLabSize),
              axis.title = element_text(size=axisTitleSize),
              plot.title = element_text(size=tabTitleSize),
              legend.title = element_text(size=tabLabSize),
              legend.text = element_text(size=tabLabSize)) + 
        guides(shape = guide_legend(title=legTitle,
                                    keywidht=tabLabSize,
                                    keyheight=tabLabSize))
    ## remove legend
    if(!legend) g2 <- g2 + theme(legend.position = "none")
    res1 <- list("table"=g2,
                "plot"=g1)
    class(res1) <- c("tableAndPlot", "list")
    return(res1)
}
## declare variables (for R CMD check)
## st1 is vector for strata identification
surv <- n.risk <- n.censor <- n.event <- upper <- lower <- NULL
.SD <- st1 <- stNames <- st <- s1 <- minT <- l <- maxT <- u <- NULL
