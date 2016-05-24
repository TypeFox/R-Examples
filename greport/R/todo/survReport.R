#' Survival Report
#'
#' Generate a survival report with Kaplan-Meier estimates,
#'
#' If both \code{group} and \code{treat.group} are present, plots will be
#' generated for each level of group.
#'
#' @param etime character. Name of column with event time data.
#' @param event character. Name of column with event status data.
#' @param treat character. Name of column with treatment data.
#' @param group character. Name of column with group data.
#' @param treat.group character. Name of column with treatment group data.
#' @param data data.frame. Data used for report.
#' @param ylabel character. Passed to \code{\link[rms]{survplot.survfit}} as
#' the \code{ylab} argument.
#' @param conf character. See \code{\link[rms]{survplot.survfit}}.
#' @param n numeric. See \code{\link[Hmisc]{ldBands}}.
#' @param labels character vector.  See \code{\link[Hmisc]{plot.ldBands}}.
#' @param previousActual Passed as the second component to the
#' \code{actual} argument in the \code{\link[Hmisc]{plot.ldBands}} function.
#' @param h numeric. Height of plot. Default is 4in. See \code{\link[Hmisc]{setps}}.
#' @param append logical. If \sQuote{TRUE} output will be appended instead of overwritten.
#' @param fileName character. File name suffix.
#' @param descrip character. Used in caption to describe the predictor of interest.
#' The default is \sQuote{treatment}.
#' @param ndescrip character. Used in longcaption to provide further detail.
#' @param \dots additional arguments, passed to \code{\link[rms]{survplot}}
#' and \code{\link[Hmisc]{ldBands}}.
#' @export
#' @examples
#' \dontrun{
#'   set.seed(47)
#'   mydata <- data.frame(time=sapply(round(rnorm(1000, 36, 5)), min, 36), treat=rep(0:1, each=500))
#'   mydata$event <- ifelse(mydata$time < 36, 1, 0)
#'   survReport("time", "event", "treat", data=mydata)
#' }

survReport <- function(etime, event, treat, group=NULL, treat.group=NULL, data,
                       ylabel='Survival Probability',
                       conf=c('bars','bands','none'),
                       n=NULL, labels=NULL, previousActual=NULL, h=4, 
                       append=FALSE, fileName="trt", descrip="treatment", ndescrip = "", ...) {
  plotName <- paste('surv-km', fileName, sep=".")
  conf <- match.arg(conf)

  if(length(group) & length(treat.group)) {
    # Base K-M plot on treat.group
    lty = rep(c(1,3), 2); col = rep(gray(c(0,.7)), each = 2); lwd = 1
    S <- Surv(data[[etime]], data[[event]])
    treat.group <- as.factor(data[[treat.group]])
    startPlot(plotName, h=h)
    survplot.survfit(survfit.formula(S ~ treat.group), n.risk=TRUE, conf=conf, lwd=lwd,
       lty=lty, col=col, ylab=ylabel, ...)
    endPlot()
    # Calculate the corresponding z for each level of group 
    for(i in levels(as.factor(data[[group]]))) {
      assign(paste("S", i, sep="."), Surv(data[data[[group]]==i, etime], data[data[[group]]==i, event]))
      assign(paste("treat", i, sep="."), data[data[[group]]==i, treat])
      j <- is.na(get(paste("S", i, sep="."))) | is.na(get(paste("treat", i, sep="."))) 
      if(any(j)) {
         assign(paste("S", i, sep="."), get(paste("S", i, sep="."))[!j,])
         assign(paste("treat", i, sep="."), get(paste("treat", i, sep="."))[!j])
      }
      assign(paste("z", i, sep="."), round(sqrt(logrank(get(paste("S", i, sep=".")), as.integer(get(paste("treat", i, sep="."))))), 2))
    }
    putFig(panel = 'surv', name = plotName,
      caption = paste('Kaplan-Meier estimates by ',  descrip, '.', sep=''), 
      longcaption = paste('Kaplan-Meier cumulative event-free probability estimates by ', paste(descrip, ndescrip, sep = " "), '.', 
        '  Based on these data, the unadjusted Cox-logrank $z$ values for ', levels(as.factor(data[[group]]))[1], ' and ', levels(as.factor(data[[group]]))[2],
        ' across treatment were calculated to be ', get(paste("z", levels(as.factor(data[[group]]))[1], sep=".")), ' and ', 
        get(paste("z", levels(as.factor(data[[group]]))[2], sep=".")), 
        ', respectively.  ',
        # Inset treatment key        
        levels(treat.group)[1], ' = solid black; ', levels(treat.group)[2], ' = dotted black; ', levels(treat.group)[3], ' = solid gray; ', 
        levels(treat.group)[4],' = dotted gray.', sep=""), append = append)
  }
  else {
    # Base everything on treat
    lwd = c(1,2); lty = c(1,1); col = gray(c(0,.7))
    S <- Surv(data[[etime]], data[[event]])
    treat <- as.factor(data[[treat]])
    if(attributes(survfit.formula(S ~ treat)$strata)$names[1]=="B") {
       col=gray(c(0.7, 0))
    }
    startPlot(plotName, h=h)
    survplot.survfit(survfit.formula(S ~ treat), n.risk=TRUE, conf=conf, lwd=lwd,
       lty=lty, col=col, ylab=ylabel, ...)
    endPlot()
    # Calculate the corresponding z 
    i <- is.na(S) | is.na(treat) 
    if(any(i)) {
      S <- S[!i,]
      treat <- treat[!i]
    }
    z <- round(sqrt(logrank(S, as.integer(treat))), 2)
    putFig(panel = 'surv', name = plotName,
      caption = paste('Kaplan-Meier estimates by ',  descrip, '.', sep=''), 
      longcaption = paste('Kaplan-Meier cumulative event-free probability estimates by ', paste(descrip, ndescrip, sep = " "), ".", 
           '  Based on these data, an unadjusted Cox-logrank $z$ value was calculated to be ', z, ".  ", 
           # Inset treatment key        
           levels(treat)[1], ':\\rule[.05in]{.25in}{.5pt}; ', levels(treat)[2], ':\\textcolor[gray]{0.7}{\\rule[.05in]{.25in}{1.25pt}}.', sep=""), 
           append = append)
  }

# LEAVE AS IS IN survReport.s
  if(length(n)) {
    p <- ldBands(n=n, ...)
    i <- is.na(S) | is.na(treat)
    if(any(is.na(i))) {
      S <- S[!i,]
      treat <- treat[!i]
    }
    z <- logrank(S, as.integer(treat))
    startPlot('surv-monitor', h=h)
    plot(p, labels=labels, actual=c(z, previousActual))
    endPlot()
    putFig('surv','surv-monitor',
           'Group-sequential monitoring boundaries',
           paste('Lan-DeMets group-sequential monitoring boundaries using the Obrien-Fleming alpha-spending function with',n,'looks equally spaced in time.  Points indicate observed Cox-logrank Z statistics.'), append=append)
  }
}
