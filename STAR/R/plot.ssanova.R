##############################################################################
##############################################################################
##############################################################################
plot.ssanova <- function(x,
                         y,
                         include,
                         ask=FALSE,
                         ncol=2,
                         nrow=3,
                         ...
                         ) {
#######################################################################
### Method / function plot.ssanova
### A plot method for ssanova objects tailored to their use in STAR
### ----------------------------------------------------------
### Arguments:
###  x: a "ssanova" object.
###  y: not used, only included for compatibility with generic method.
###  include: a character string with the model terms one wants to plot.
###           if missing all terms are plotted. 
###  ask: a logical. If TRUE terms are plotted (on a common y scale) one
###       one after the other and the user is invited to hit the "enter"
###       key to generate the next plot. If FALSE (default) all terms
###       are drawn on a suitable number of X11 devices. The number of
###       terms on each device is controlled by arguments ncol and nrow.
###  ncol: the number of columns of the display matrix used on each device
###        when "ask" is set to FALSE.
###  nrow: the number of lines of the display matrix used on each device
###        when "ask" is set to FALSE.
### ...:  not used only there for method definition compatibility.
#######################################################################
  ## Identify terms
  allTerms <- x$terms$labels[-1]
  ## Keep only terms made of at most interactions between 2 variables
  allTerms <- allTerms[sapply(strsplit(allTerms,":"),length)<3]

  if (missing(include)) {
    include <- allTerms
  } else {
    include <- include[include %in% allTerms]
    if (length(include) == 0) stop("include elements incorrect.")
  }

  nbPlots <- length(include)
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  } else {
    nbPanels <- (nbPlots-1) %/% (ncol*nrow) + 1
    panelIdx <- 1
    layout(matrix(1:(ncol*nrow),nrow=nrow,ncol=ncol))
  }

  allPred <- lapply(seq(along=include),
                    function(plotIdx) {
                      theTerm <- include[plotIdx]
                      isInteraction <- length(strsplit(theTerm,":")[[1]]) == 2
                      if (!isInteraction) {
                        quickPredict(x,include=theTerm,length.out=501)
                      } else {
                        quickPredict(x,include=theTerm,length.out=101)
                      }
                    }
                    )

  yMin <- min(sapply(allPred, function(l) min(l$est.mean)))
  yMax <- max(sapply(allPred, function(l) max(l$est.mean)))
  yDelta <- yMax-yMin
  ylim <- c(yMin-0.02*yDelta,yMax+0.02*yDelta)
  
  for (plotIdx in seq(along=include)) {
    theTerm <- include[plotIdx]
    isInteraction <- length(strsplit(theTerm,":")[[1]]) == 2
    if (!isInteraction) plot(allPred[[plotIdx]],
                             ylim=ylim,
                             panel.first=grid(col=1))
    else {
      myLevels <- unique(sort(c(0,seq(ylim[1],ylim[2],len=10))))
      myLabels <- round(myLevels,digits=2)
      contour(allPred[[plotIdx]],what="mean",
              levels=myLevels,col=2,lwd=2,labcex=1,labels=myLabels)
      contour(allPred[[plotIdx]],what="sd",
              nlevels=5,col=1,lty=3,add=TRUE,labcex=1)
    }
    if (!ask) {
      if (plotIdx %% (ncol*nrow) == 0 && panelIdx < nbPanels) {
        X11()
        layout(matrix(1:(ncol*nrow),nrow=nrow,ncol=ncol))
        panelIdx <- panelIdx + 1
      }
    } ## End of conditional on !ask
  } ## End of loop on plotIdx
  
}
##############################################################################
##############################################################################
##############################################################################
plot.ssanova0 <- plot.ssanova
