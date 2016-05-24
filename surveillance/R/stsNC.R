######################################################################
# initialize-method for "stsNC" objects
######################################################################

init.stsNC <- function(.Object, ...,
                       reportingTriangle, predPMF, pi, truth, delayCDF, SR)
{
    .Object <- callNextMethod()  # use initialize,sts-method

    ## initialize defaults for extra stsNC-slots or check supplied values
    dimObserved <- dim(.Object@observed)
    if (missing(pi)) {
        .Object@pi <- array(NA_integer_, dim = c(dimObserved, 2L))
    } else {
        dimPI <- dim(.Object@pi)
        if (length(dimPI) != 3 || any(dimPI != c(dimObserved, 2L)))
            stop("dim(pi) = (", paste0(dimPI, collapse=","), ")")
    }
    if (missing(SR)) {
        .Object@SR <- array(NA_real_, dim = c(nrow(.Object@observed),0L,0L))
    } else {
        stopifnot(length(dim(.Object@SR)) == 3)
    }
    if (missing(truth))
        .Object@truth <- as(.Object, "sts")
  
    return(.Object)
}

setMethod("initialize", "stsNC", init.stsNC)


######################################################################
# Special coerce method to account for consistent dimensions
######################################################################

setAs(from = "sts", to = "stsNC", function (from) {
    new("stsNC", from,
        pi = array(NA_real_, dim = c(dim(from@observed), 2L)),
        truth = from,
        SR = array(NA_real_, dim = c(nrow(from@observed), 0L, 0L)))
})


######################################################################
# plot-method for the "stsNC" class, which starts by
# using the inherited method, but with some additional plotting
# put into the .hookFunSpecial function.
#
# Parameters:
#  same as the for the plot method of sts objects.
######################################################################

setMethod(f="plot", signature=signature(x="stsNC", y="missing"),
          function (x, type = observed ~ time | unit, ...) {

              ## environment of hook function will be set to evaluation 
              ## environment of stsplot_time1() and only then be called
              legend.opts <- lty <- lwd <-
                  "accommodate tools:::.check_code_usage_in_package()"
              
              #Hook function specifically for nowcasting objects.
              nowcastPlotHook <- function() {
                  #Define some colors for the plotting as well as some plot symbols
                  color <- surveillance.options("colors")
                  pchList   <- c(nowSymbol=10)
                  
                  #Prolong line of last observation (this should go into the plot function
                  idx <- nrow(x) - which.max(!is.na(rev(upperbound(x)))) + 1
                  #Continue line from plot - use same style as stsplot_time1
                  lines( idx+c(-0.5,0.5), rep(upperbound(x)[idx,],2),col=col[3],lwd=lwd[3],lty=lty[3])

                  #Add the prediction intervals as bars (where not NA). Conf level
                  #is found in x@control$alpha
                  idxt <- which(apply(x@pi[1:nrow(x),1,],1,function(x) all(!is.na(x))))
                  for (i in idxt) {
                      lines( i+c(-0.3,0.3), rep(x@pi[i,,1],2),lty=1,col=color["piBars"])
                      lines( i+c(-0.3,0.3), rep(x@pi[i,,2],2),lty=1,col=color["piBars"])
                      lines( rep(i,each=2), x@pi[i,,],lty=2,col=color["piBars"])
                  }

                  #Extract now date and date range of the plotting
                  startDate <- epoch(x)[1]

                  #Add "now" symbol on x-axis. Plotting now takes possible temporal aggregation into account.
                  #points(x@control$now-startDate+1,0,pch=pchList["nowSymbol"],col=color["nowSymbol"],cex=1.5)
                  points(x@control$timeDelay(startDate,x@control$now)+1,0,pch=pchList["nowSymbol"],col=color["nowSymbol"],cex=1.5)
                  #Add this to the legend
                  if (!is.null(legend.opts)) {
                      legend(x="topright",c("Now"),pch=pchList["nowSymbol"],col=color["nowSymbol"],bg="white")
                  }
                  
                  return(invisible())
              }

              callNextMethod(x=x, type=type, ..., .hookFuncInheritance=nowcastPlotHook)
          })



