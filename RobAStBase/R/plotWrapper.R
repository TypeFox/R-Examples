##########################################
##                                      ##
##    Wrapper for infoPlot.R            ##
##    (infoPlot method for IC)          ##
##                                      ##
##########################################

##############################################################
#' Merging Lists
#'
#' \code{.merge.lists} takes two lists and merges them.
#'
#' @param a the first list
#'
#' @param b the second list
#'
#' @return the merged list
#'
#' @keywords internal
#' @rdname mergelists
#'
##############################################################

### aditional function
.merge.lists <- function(a, b){
  a.names <- names(a)
  b.names <- names(b)
  m.names <- sort(unique(c(a.names, b.names), fromLast = TRUE))
  sapply(m.names, function(i) {
    if (is.list(a[[i]]) & is.list(b[[i]])) .merge.lists(a[[i]], b[[i]])
    else if (i %in% b.names) b[[i]]
    else a[[i]]
  }, simplify = FALSE)
}

##############################################################
#' Wrapper function for information plot method
#'
#' The wrapper takes most of arguments to the plot method
#' by default and gives a user possibility to run the
#' function with low number of arguments
#'
#' @param IC object of class \code{IC}
#'
#' @param data optional data argument --- for plotting observations into the plot
#'
#' @param ... additional parameters (in particular to be passed on to \code{plot})
#'
#' @param alpha.trsp the transparency argument (0 to 100) for ploting the data
#'
#' @param with.legend the flag for showing the legend of the plot
#'
#' @param rescale the flag for rescaling the axes for better view of the plot
#'
#' @param withCall the flag for the call output
#'
#' @return invisible(NULL)
#
#' @section Details: Calls \code{infoPlot} with suitably chosen defaults. If \code{withCall == TRUE}, the call to \code{infoPlot} is returned
#'
#' @export
#' @rdname InfoPlotWrapper
#'
#'
#' @examples
#' # Gamma
#' fam  <-  GammaFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y <- distribution(fam)
#' data  <-  r(Y)(1000)
#' InfoPlot(IC, data, withCall = FALSE)
#'
##############################################################

##IC - influence curve
##data - dataset
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
InfoPlot <- function(IC, data,...,alpha.trsp = 100,with.legend = TRUE, rescale = FALSE ,withCall = TRUE){
  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  if(missing(IC)) stop("Argument 'IC' must be given as argument to 'InfoPlot'")
  if(missing(data)) data <- NULL
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(missing(data)){
    alpha.trsp <- 100
  } else {
    if(is.null(mc$alpha.trsp)){
      alpha.trsp <- 30
      if(length(data) < 1000){
        alpha.trsp <- 50
      }
      if(length(data) < 100){
        alpha.trsp <- 100
      }
    }
  }
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$rescale)) mc$rescale <- FALSE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##

  ## Scaling of the axes
  scaleList <- rescaleFunction(eval(IC@CallL2Fam), FALSE, mc$rescale)

  argsList <- c(list(object = substitute(IC)
                     ,data = substitute(data)
                     ,withSweave = substitute(getdistrOption("withSweave"))
                     ,lwd = substitute(par("lwd"))
                     ,lty = substitute("solid")
                     ,colI = substitute(grey(0.5))
                     ,lwdI = substitute(0.7*par("lwd"))
                     ,ltyI = substitute("dotted")
                     ,main = substitute(FALSE)
                     ,inner = substitute(TRUE)
                     ,sub = substitute(FALSE)
                     ,col.inner = substitute(par("col.main"))
                     ,cex.inner = substitute(0.8)
                     ,bmar = substitute(par("mar")[1])
                     ,tmar = substitute(par("mar")[3])
                     ,with.legend = substitute(TRUE)
                     ,legend = substitute(NULL)
                     ,legend.bg = substitute("white")
                     ,legend.location = substitute("bottomright")
                     ,legend.cex = substitute(0.8)
                     ,scaleN = substitute(9)
                     ,mfColRow = substitute(TRUE)
                     ,to.draw.arg = substitute(NULL)
                     ,cex.pts = substitute(1)
                     ,col.pts = substitute(addAlphTrsp2col(rgb(0,255,0,maxColorValue=255), substitute(alpha.trsp)))
                     ,pch.pts = substitute(19)
                     ,jitter.fac = substitute(1)
                     ,with.lab = substitute(FALSE)
                     ,lab.pts = substitute(NULL)
                     ,lab.font = substitute(NULL)
                     ,alpha.trsp = substitute(NA)
                     ,which.lbs = substitute(NULL)
                     ,which.Order  = substitute(NULL)
                     ,return.Order = substitute(FALSE)
                     ,ylab.abs = substitute("absolute information")
                     ,ylab.rel= substitute("relative information")
                     ,adj = substitute(0.5)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("o")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue")
    ), scaleList)

  ##parameter for plotting
  if(mc$with.legend)
  {
    argsList$col.main <- "black"
    argsList$col.lab <- "black"
  }
  else
  {
    argsList$col.main <- "white"
    argsList$col.lab <- "white"
  }

  args <- .merge.lists(argsList, dots)
  ###
  ### 3. build up the call but grab it and write it into an object
  ###
  cl <- substitute(do.call(infoPlot,args0), list(args0=args))
  ### manipulate it so that the wrapper do.call is ommitted
  cl0 <- as.list(cl)[-1]
  mycall <- c(cl0[1],unlist(cl0[-1]))
  mycall <- as.call(mycall)
  ###
  ### 4. evaluate the call (i.e., produce the graphic)
  ###
  eval(mycall)
  ###
  ### 5. return the call (if withCall==TRUE)
  ###
  if(mc$withCall) print(mycall)

}


################################################################################

##########################################
##                                      ##
##    Wrapper for AllPlot.R             ##
##    (plot method for IC)              ##
##                                      ##
##########################################


##############################################################
#' Wrapper function for plot method for IC
#'
#' The wrapper takes most of arguments to the plot method
#' by default and gives a user possibility to run the
#' function with low number of arguments
#'
#' @param IC object of class \code{IC}
#'
#' @param y optional data argument --- for plotting observations into the plot
#'
#' @param ... additional parameters (in particular to be passed on to \code{plot})
#'
#' @param alpha.trsp the transparency argument (0 to 100) for ploting the data
#'
#' @param with.legend the flag for showing the legend of the plot
#'
#' @param rescale the flag for rescaling the axes for better view of the plot
#'
#' @param withCall the flag for the call output
#'
#' @return invisible(NULL)
#
#' @section Details: Calls \code{plot} with suitably chosen defaults; if \code{withCall == TRUE}, the call to \code{plot} is returned.
#'
#' @export
#' @rdname PlotICWrapper
#'
#' @examples
#' # Gamma
#' fam <- GammaFamily()
#' rfam <- InfRobModel(fam, ContNeighborhood(0.5))
#' IC <- optIC(model = fam, risk = asCov())
#' Y <- distribution(fam)
#' y <- r(Y)(1000)
#' PlotIC(IC, y, withCall = FALSE)
##############################################################

##IC - influence curve
##y - dataset
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
PlotIC <- function(IC, y,...,alpha.trsp = 100, with.legend = TRUE, rescale = FALSE ,withCall = TRUE){
  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  if(missing(IC)) stop("Argument 'IC' must be given as argument to 'PlotIC'")
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(missing(y)){
    alpha.trsp <- 100
  } else {
    if(is.null(mc$alpha.trsp)){
      alpha.trsp <- 30
      if(length(y) < 1000){
        alpha.trsp <- 50
      }
      if(length(y) < 100){
        alpha.trsp <- 100
      }
    }
  }
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$rescale)) mc$rescale <- FALSE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##

  ## Scaling of the axes
  scaleList <- rescaleFunction(eval(IC@CallL2Fam), !missing(y), mc$rescale)

  argsList <- c(list(x = substitute(IC)
                     ,withSweave = substitute(getdistrOption("withSweave"))
                     ,main = substitute(FALSE)
                     ,inner = substitute(TRUE)
                     ,sub = substitute(FALSE)
                     ,col.inner = substitute(par("col.main"))
                     ,cex.inner = substitute(0.8)
                     ,bmar = substitute(par("mar")[1])
                     ,tmar = substitute(par("mar")[3])
                     ,with.legend = substitute(FALSE)
                     ,legend = substitute(NULL)
                     ,legend.bg = substitute("white")
                     ,legend.location = substitute("bottomright")
                     ,legend.cex = substitute(0.8)
                     ,withMBR = substitute(FALSE)
                     ,MBRB = substitute(NA)
                     ,MBR.fac = substitute(2)
                     ,col.MBR = substitute(par("col"))
                     ,lty.MBR = substitute("dashed")
                     ,lwd.MBR = substitute(0.8)
                     ,scaleN = substitute(9)
                     ,mfColRow = substitute(TRUE)
                     ,to.draw.arg = substitute(NULL)
                     ,adj = substitute(0.5)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("o")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue")
    ), scaleList)
  if(missing(y)){c(argsList, y = substitute(y)
                     ,cex.pts = substitute(0.3)
                     ,col.pts = substitute(addAlphTrsp2col(rgb(0,255,0,maxColorValue=255), substitute(alpha.trsp)))
                     ,pch.pts = substitute(19)
                     ,jitter.fac = substitute(1)
                     ,with.lab = substitute(FALSE)
                     ,lab.pts = substitute(NULL)
                     ,lab.font = substitute(NULL)
                     ,alpha.trsp = substitute(NA)
                     ,which.lbs = substitute(NULL)
                     ,which.Order  = substitute(NULL)
                     ,return.Order = substitute(FALSE)
                     ,scaleN = substitute(9)
                     ,adj = substitute(0.5)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("o")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue"))
  }

  ##parameter for plotting
  if(mc$with.legend)
  {
    argsList$col.main <- "black"
    argsList$col.lab <- "black"
  }
  else
  {
    argsList$col.main <- "white"
    argsList$col.lab <- "white"
  }

  args <- .merge.lists(argsList, dots)
  ###
  ### 3. build up the call but grab it and write it into an object
  ###
  cl <- substitute(do.call(plot,args0), list(args0=args))
  ### manipulate it so that the wrapper do.call is ommitted
  cl0 <- as.list(cl)[-1]
  mycall <- c(cl0[1],unlist(cl0[-1]))
  mycall <- as.call(mycall)
  ###
  ### 4. evaluate the call (i.e., produce the graphic)
  ###
  eval(mycall)
  ###
  ### 5. return the call (if withCall==TRUE)
  ###
  if(mc$withCall) print(mycall)

}

################################################################################

##########################################
##                                      ##
##    Wrapper for comparePlot)          ##
##                                      ##
##########################################


##############################################################
#' Wrapper function for function comparePlot
#'
#' The wrapper takes most of arguments to function comparePlot
#' by default and gives a user possibility to run the
#' function with low number of arguments
#'
#' @param IC1 object of class \code{IC}
#'
#' @param IC2 object of class \code{IC}
#'
#' @param IC3 object of class \code{IC}
#'
#' @param IC4 object of class \code{IC}
#'
#' @param y optional data argument --- for plotting observations into the plot
#'
#' @param ... additional parameters (in particular to be passed on to \code{plot})
#'
#' @param alpha.trsp the transparency argument (0 to 100) for ploting the data
#'
#' @param with.legend the flag for showing the legend of the plot
#'
#' @param rescale the flag for rescaling the axes for better view of the plot
#'
#' @param withCall the flag for the call output
#'
#' @return invisible(NULL)
#
#' @section Details: Calls \code{comparePlot} with suitably chosen defaults; if \code{withCall == TRUE}, the call to \code{comparePlot} is returned.
#'
#'
#' @export
#' @rdname ComparePlotWrapper
#'
#' @examples
#' # Gamma
#' fam <- GammaFamily()
#' rfam <- InfRobModel(fam, ContNeighborhood(0.5))
#' IC1 <- optIC(model = fam, risk = asCov())
#' IC2 <- makeIC(list(function(x)sin(x),function(x)x^2), L2Fam = fam)
#' Y <- distribution(fam)
#' y <- r(Y)(1000)
#' ComparePlot(IC1, IC2, y, withCall = TRUE)
##############################################################

##IC - influence curve
##y - dataset
## with.legend - optional legend indicator
## withCall - optional indicator of the function call
#
ComparePlot <- function(IC1, IC2, y, ..., IC3=NULL, IC4=NULL,
        alpha.trsp = 100, with.legend = TRUE, rescale = FALSE ,withCall = TRUE){

  ###
  ### 1. grab the dots (and manipulate it within the wrapper function)
  ###
  ###
  ### do something to fix the good default arguments
  ###
  if(missing(IC1)) stop("Argument 'IC1' must be given as argument to 'ComparePlot'")
  if(missing(IC2)) stop("Argument 'IC2' must be given as argument to 'ComparePlot'")
  mc <- as.list(match.call(expand.dots = FALSE))[-1]
  dots <- mc$"..."
  if(missing(y)){
    alpha.trsp <- 100
  } else {
    if(is.null(mc$alpha.trsp)){
      alpha.trsp <- 30
      if(length(y) < 1000){
        alpha.trsp <- 50
      }
      if(length(y) < 100){
        alpha.trsp <- 100
      }
    }
  }
  if(is.null(mc$with.legend)) mc$with.legend <- TRUE
  if(is.null(mc$rescale)) mc$rescale <- FALSE
  if(is.null(mc$withCall)) mc$withCall <- TRUE
  iny <- if(missing(y)) TRUE else is.null(y)
  ###
  ### 2. build up the argument list for the (powerful/fullfledged)
  ### graphics/diagnostics function;
  ##

  ## Scaling of the axes
  scaleList <- rescaleFunction(eval(IC1@CallL2Fam), iny, rescale)

  argsList <- .merge.lists(list(obj1 = substitute(IC1)
                     ,obj2 = substitute(IC2)
                     ,obj3 = NULL
                     ,obj4 = NULL
                     ,forceSameModel = FALSE
                     ,data = NULL
                     ,withSweave = substitute(getdistrOption("withSweave"))
                     ,main = substitute(FALSE)
                     ,inner = substitute(TRUE)
                     ,sub = substitute(FALSE)
                     ,col.inner = substitute(par("col.main"))
                     ,cex.inner = substitute(0.8)
                     ,bmar = substitute(par("mar")[1])
                     ,tmar = substitute(par("mar")[3])
                     ,with.legend = substitute(FALSE)
                     ,legend = substitute(NULL)
                     ,legend.bg = substitute("white")
                     ,legend.location = substitute("bottomright")
                     ,legend.cex = substitute(0.8)
                     ,withMBR = substitute(FALSE)
                     ,MBRB = substitute(NA)
                     ,MBR.fac = substitute(2)
                     ,col.MBR = substitute(par("col"))
                     ,lty.MBR = substitute("dashed")
                     ,lwd.MBR = substitute(0.8)
                     ,scaleX = FALSE
                     ,scaleX.fct = NULL
                     ,scaleX.inv = NULL
                     ,scaleY = FALSE
                     ,scaleY.fct = pnorm
                     ,scaleY.inv=qnorm
                     ,scaleN = 9
                     ,x.ticks = NULL
                     ,y.ticks = NULL
                     ,mfColRow = substitute(TRUE)
                     ,to.draw.arg = substitute(NULL)
                     ,cex.pts = substitute(1)
                     ,col.pts = substitute(c(1,2,3,4))
                     ,pch.pts = substitute(19)
                     ,jitter.fac = substitute(1)
                     ,with.lab = substitute(FALSE)
                     ,lab.pts = substitute(NULL)
                     ,lab.font = substitute(NULL)
                     ,alpha.trsp = substitute(alpha.trsp)
                     ,which.lbs = substitute(NULL)
                     ,which.Order  = substitute(NULL)
                     ,return.Order = substitute(FALSE)
                     ,adj = substitute(0.5)
                     ,cex.main = substitute(1.5)
                     ,cex.lab = substitute(1.5)
                     ,cex = substitute(1.5)
                     ,bty = substitute("o")
                     ,panel.first= substitute(grid())
                     ,col = substitute("blue")
    ), scaleList)
    
    if(!is.null(IC3)) argsList$obj3 <- substitute(IC3)
    if(!is.null(IC4)) argsList$obj4 <- substitute(IC4)

    if(!missing(y))  argsList$data <- substitute(y)

  ##parameter for plotting
  if(mc$with.legend)
  {
    argsList$col.main <- "black"
    argsList$col.lab <- "black"
  }
  else
  {
    argsList$col.main <- "white"
    argsList$col.lab <- "white"
  }

  args <- .merge.lists(argsList, dots)
  wn <- which(names(args) %in% c("obj1", "obj2"))
  args <- c(args[wn],args[-wn])
  ###
  ### 3. build up the call but grab it and write it into an object
  ###
  cl <- substitute(do.call(comparePlot,args0), list(args0=args))
  ### manipulate it so that the wrapper do.call is ommitted
  cl0 <- as.list(cl)[-1]
  mycall <- c(cl0[1],unlist(cl0[-1]))
  mycall <- as.call(mycall)
  ###
  ### 4. evaluate the call (i.e., produce the graphic)
  ###
  eval(mycall)
  ###
  ### 5. return the call (if withCall==TRUE)
  ###
  if(mc$withCall) print(mycall)

}

