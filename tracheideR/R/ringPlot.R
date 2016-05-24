#' @name ringPlot
#' @title Plot intra-ring variation of tracheid features
#' (with the possibility of plotting a climatic variable)
#' @description
#' This function plots the intra-ring variation of a tracheid feature (e.g. LD, CWT or LD/CWT)
#' along the growing season and, optionally, the intra-annual variation of a climatic variable can be added to the same plot.
#' @param traq ordered sequence of a tracheid feature (e.g. LD)
#' @param varMean vector with monthly values of a given environmental variable
#' @param varYear vector with monlthy values of an environmental variable for a specific year
#' @param m0 the start of the growing season (in months); the default value is 3 that corresponds to day of the year (doy) 60.
#' @param mt the moment of transition from earlywood to latewood; the default value is 6.75.
#' @param m1 the end of the growing season; default value is 11.
#' @param type a string that defines the tracheid features to be plotted, defaults "LD"
#' @param ylab the y axis title (for the tracheid feature variable), default value is \code{type}
#' @param main an overall title for the plot, if no \code{string} is supplied no title is added to the plot
#' @param addGS logical; if \code{TRUE} add the growing season length to the axis 1, defaults \code{TRUE}
#' @param addMonths logical; if \code{TRUE} add months to the axis 1, defaults \code{TRUE}
#' @param ... graphical parameters for \code{plot} may also be passed as arguments to this function
#' @param varMeanCol the default value, \code{"grey80"}, gives the color to plot the mean environmental variable
#' @param varYearCol the color to plot monthly environmental values; the default value is \code{"red"}
#' @param varName the y axis title (for the environmental variable),default value is \code{""}; if no \code{string} is supplied no title is added to environmental axis
#' @details This function returns an invisible \code{data.frame} (used to produced the graph)
#' @export
#' @importFrom tgram tgram
#' @importFrom graphics axis lines mtext par plot title
#' @importFrom grDevices extendrange
#' @examples
#'
#' ## Not run:
#'
#' # year 2010
#' y2010raw <- getTrac(tch$y2010, scale = .169)
#' y2010std <- tracheider(y2010raw)
#' par(oma = c(2,1,1,0.5))
#' par(mar = c(2,4,1,4))
#' y2010LD <- ringPlot(traq = y2010std, varMean = colMeans(sw),
#' varYear = sw["2010",], main=2010,type = "LD",  ylim = c(0,45),
#' ylab = expression(paste("LD (", mu,"m)")),varName = "Soil moisture")
#'
#' # year 2013
#' y2013raw <- getTrac(tch$y2013, scale = .169)
#' y2013std <- tracheider(y2013raw)
#' y2013LD <- ringPlot(traq = y2013std, varMean = colMeans(sw),
#' varYear = sw["2013",], main=2013,type = "LD",  ylim = c(0,45),
#' ylab = expression(paste("LD (", mu,"m)")),varName = "Soil moisture")
#'
#' # 2010 & 2013 in the same plot
#' par(mfcol = c(2,1))
#' par(oma = c(2,1,1,0.5))
#' par(mar = c(2,4,1,4))
#'
#' ringPlot(y2010std,  varMean = colMeans(sw), varYear = sw["2010",],
#' type = "LD", ylab = "", main=2010, addGS = FALSE, addMonths = FALSE)
#'
#' ringPlot(y2013std,  varMean = colMeans(sw), varYear = sw["2013",],
#' type = "LD", ylab=expression(paste("LD (", mu,"m)")),
#' main = 2013, addGS = TRUE, varName= "Soil moisture")
#'
#' ## End(not run)


ringPlot <-
  function (traq, varMean = NULL, varYear = NULL, m0 = 3, mt = 6.75, m1 = 11,
            type = c("LD","CWT","LWratio"),
            ylab = match.arg(type), main = "", addGS = TRUE, addMonths = addGS, ...,
            varMeanCol ="grey80", varYearCol = "red", varName = "") {



    type <- match.arg(type)


    traq$LD <- rowMeans(data.frame(traq$LD))
    traq$CWT <- rowMeans(data.frame(traq$CWT))
    traq$LWratio <- rowMeans(data.frame(traq$LWratio))

    nE <- min(which(traq$LD / traq$CWT <= 2))-1

    if (type == "LWratio") {
      nE <- min(which(traq$LWratio <= 2))-1
    }

    day2month = function (day) {
      md <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28)
      mc <- cumsum(md)
      mc0 <- c(0,mc)
      m <- 0
      while ((day - mc[m + 1]) >= 0) {
        m = m + 1
      }
      m + (day - mc0[m + 1]) / md[m + 1]

    }

    if (max(m1) > 12) {
      m0 <- day2month(m0)
      mt <- day2month(mt)
      m1 <- day2month(m1)
    }
    y <- switch(type,
                "LD" = traq$LD,
                "CWT" = traq$CWT,
                "LWratio" = traq$LWratio)


    f2t <- function (y, E.start, E.end = NULL, L.end,  nE = NULL) {
      n = length(y)
      if (is.null(E.end)) {
        n1 = (L.end - E.start) / (n - 1)
        x = E.start + 0:(n - 1) * n1
        pch = rep(1,n)
        if (!is.null(nE))
          pch = c(rep(1,nE), rep(16,n - nE))
        return(data.frame(x,y,pch))
      }
      nL = n - nE
      n1 = (E.end - E.start) / (nE - 1)
      x1 = E.start + 0:(nE - 1) * n1
      n2 = (L.end - E.end - n1) / (nL - 1)
      x2 = E.end + n1 + 0:(nL - 1) * n2
      x = c(x1,x2)
      pch = c(rep(1,nE), rep(16,nL))
      return(data.frame(x,y, pch))
    }

    PLOT <- function (x, ...) plot(x=x[,1]-0.5,y=x[,2], pch=x[,3],...)

    f2t(y,m0,mt,m1,nE)->out
    PLOT(out, ann=FALSE, axes=FALSE,xlim=c(0.5,13.5)-0.5,  xaxs="i",...)
    title(main=main)
    axis(side=1, at=1:13-.5, labels = FALSE,tcl=-0.5)


    if (addGS){
      axis(side = 1, at=c(m0,m1)-0.5, line=2, labels=FALSE, tcl=-0.6)
      mtext(side=1, line=2.3,at=6.5, "Growing season")
      axis(side = 1, at=c(m0,m1)-0.5, line=3.8, labels=FALSE, tcl=0.6)
    }

    if (addMonths){
      mtext(side = 1, line=0.5,at =seq(1.5,12.5,1)-.5, strtrim(month.abb,1))
    }

    axis(side = 2)
    mtext(side = 2, text = ylab, outer = TRUE, line = -1)
    mtext(side = 4, text = varName, outer = TRUE, line = -1)


    if (!is.null(varMean )|!is.null(varYear)){
      thicks <-pretty(extendrange(c(varMean,varYear), f=.1))
      ylim <- range(thicks)
      par(new = TRUE)
      plot(varMean, type = "l", ann = FALSE,axes = FALSE,
           xlim = c(0,13),xaxs = "i", ylim = ylim, yaxs = "i", col = varMeanCol)
      lines(x = 1:12, y = varYear, col = varYearCol)
      axis(side = 4, at = thicks, labels = formatC(thicks))
    }

    return(invisible(out))

  }

