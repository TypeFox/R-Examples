#' plotR_hist shows the histogram of temperatures with set of nests and the R function superimpose
#' @title Shows the histogram of temperatures with set of nests and the R function superimpose
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Nothing
#' @param x Result data
#' @param ... Parameters used by hist or plotR functions
#' @param ylimH Scale of histogram using ylimH=c(min, max)
#' @param ylabH Label for histogram scale
#' @description Shows the histogram of temperatures with set of nests and the R function superimpose
#' plotR_hist(data)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p)
#' plotR_hist(resultNest_4p)
#' }
#' @export

plotR_hist <- function(x, ..., ylimH=NULL, ylabH="Frequency of temperatures") {

def.par <- par(no.readonly = TRUE) # save default, for resetting...
  
  
if (class(x)=="list") {
	stop("Only one series can be used with plotR_hist()")
}

nids <- x

par(mar = c(def.par[["mar"]][1:3], 5.1))

L <- modifyList(list(result=x), list(...))

a <- do.call(plotR, L) 

ylim <- ScalePreviousPlot()$ylim[c("begin", "end")]

par(new=TRUE)

L <- modifyList(list(x=x), list(...))

L["ylim"] <- NULL
L["SE"] <- NULL
L["legend"] <- NULL

if (!is.null(ylimH)) L <- c(L, list(ylim=ylimH))

lp <- list("show.box", "set.par", "parameters", "fixed.parameters", "legend", "size", 
"scaleY", "replicate.CI", "ylabH", "xlimR", "ltyCI", "lwdCI")
L <- L[!(names(L) %in% lp)]

# x2 <- (par("usr")[1]+par("usr")[2]*26)/27
# x1 <- x2*26-par("usr")[2]/0.04
xlim <- ScalePreviousPlot()$xlim[c("begin", "end")]

L <- modifyList(L, list(xlab="", ylab="", main="", axes=FALSE, freq=FALSE, 
	xlim=xlim)) 

a <- do.call(hist, L) 

axis(side=4, ylim=par("yaxp")[1:2], las=1)
mtext(ylabH, side=4, line=3, cex=par("cex"))
par(new=TRUE)
# je rtablis l'chelle des y  celle de R
plot(x = 1, y=1, ylim=ylim, xlim=xlim, xlab="", ylab="", axes=FALSE, bty="n", type="n")

# par(def.par)  #- reset to default

}
