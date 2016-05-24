as.pyramidLikert <- function(x, ...,
                             panel.width=.48,
                             px=list(
                               L=c(0, panel.width),
                               R=c(1-panel.width, 1),
                               M=c(panel.width, 1-panel.width)),
                             keepLegend=(length(x$legend$bottom$args$text) > 2),
                             xlab.top=list(
                               L=list(x$legend$bottom$args$text[1]),
                               R=list(x$legend$bottom$args$text[2]),
                               M=list(x$ylab, just=1))) {
  if (!inherits(x, "trellis") || length(x$panel.args) > 1)
    stop("'as.pyramidLikert' requires a single-panel 'trellis' object.", call.=FALSE)
  class(x) <- c("pyramidLikert", class(x))
  
  if (missing(xlab.top) && !is.null(x$xlab.top)) {
    xlab.top$L <- x$xlab.top[1]
    xlab.top$R <- x$xlab.top[2]
  }
  attr(x, "xlab.top") <- xlab.top
  attr(x, "px") <- px
  attr(x, "keepLegend") <- keepLegend
  x
}

print.pyramidLikert <- function(x, ...,
                                 panel.width=.48,
                                 px=list(
                                   L=c(0, panel.width),
                                   R=c(1-panel.width, 1),
                                   M=c(panel.width, 1-panel.width)),
                                keepLegend=(length(x$legend$bottom$args$text) > 2),
                                xlab.top=list(
                                  L=list(x$legend$bottom$args$text[1]),
                                  R=list(x$legend$bottom$args$text[2]),
                                  M=list(x$ylab, just=1))) {

  if (!is.null(attr(x, "xlab.top")))
    xlab.top <- attr(x, "xlab.top")
  else {
    if (missing(xlab.top) && !is.null(x$xlab.top)) {
      xlab.top$L <- x$xlab.top[1]
      xlab.top$R <- x$xlab.top[2]
    }}

  if (missing(keepLegend) && !is.null(attr(x, "keepLegend")))
    keepLegend <- attr(x, "keepLegend")

  if (missing(panel.width) && missing(px) && !is.null(attr(x, "px")))
    px <- attr(x, "px")
    
    
    
  ## x <- plot.likert(x, ...)
  ## x is a single-panel trellis object
  if (length(x$panel.args) > 1)
    stop("'print.pyramidLikert' requires a single-panel 'trellis' object.", call.=FALSE)
  K <- x
  class(K) <- class(K)[class(K) != "pyramidLikert"]
  
  K.legend <- K$legend
  if (keepLegend)
    K$legend$bottom$args <- emptyLegend(K$legend$bottom$args)
  else
    K$legend <- NULL
  K$x.limits <- max(abs(K$x.limits)) * c(-1,1)
  K.ylab <- ifelse(is.null(K$ylab), " ", K$ylab)
  K$ylab <- NULL
  K <- update(K, par.settings=list(
                   layout.heights=list(
                     main.key.padding=2.5,
                     key.axis.padding=0,
                     axis.top=.75)))
  
  L <- K
  L$x.limits <- c(K$x.limits[1], 0)
  L$panel.args[[1]]$x[K$panel.args[[1]]$x > 0] <- 0
  L$xlab.top <- xlab.top[["L"]]
  L$main <- " "
  if (!is.null(L$sub)) L$sub <- " "
  L$y.scales$labels <- NULL
  L$xlab <- " "
  L <- update(L, par.settings=list(
                   layout.widths=list(
                     axis.right=0,       
                     axis.key.padding=0, 
                     ylab.right=0,       
                     key.right=0,        
                     right.padding=0)))
  L$x.limits <- -L$x.limits
  L$panel.args[[1]]$x <- -L$panel.args[[1]]$x
  
  R <- K
  R$x.limits <- c(0, K$x.limits[2])
  R$panel.args[[1]]$x[K$panel.args[[1]]$x < 0] <- 0
  R$xlab.top <- xlab.top[["R"]]
  R$main <- " "
  if (!is.null(R$sub)) R$sub <- " "
  R$y.scales$labels <- NULL
  R$xlab <- " "
  R  <- update(R, par.settings=list(
                    layout.widths=list(
                      left.padding=0,     
                      key.left=0,         
                      key.ylab.padding=0, 
                      ylab=0,             
                      ylab.axis.padding=0,
                      axis.left=0)))

  
  M <- K
  ##            scales=list(y=list(tck=0)),
  M$par.settings$axis.line$col <- "transparent"
  M$panel <- function(...) {}
  M$xlab.top <- xlab.top[["M"]]
  M$x.limits[] <- " "
  M <- update(M, par.settings=list(
                    layout.widths=list(
                      left.padding=1,     
                      key.left=0,         
                      key.ylab.padding=0, 
                      ylab=0,             
                      ylab.axis.padding=0,
                      axis.left=1,        
                      axis.panel=0,       
                      strip.left=0,       
                      panel=1,            
                      between=0,          
                      axis.right=0,       
                      axis.key.padding=0, 
                      ylab.right=0,       
                      key.right=0,        
                      right.padding=1
                      )))
  M$plot.args$panel.width <- list(x=0.001, units="mm", data=NULL)
  M$x.scales$labels[] <- ""
  if (keepLegend) M$legend$bottom$args <- x$legend$bottom$args

  ## if (keepLegend) {
  ##   M.width <- M
  ##   M.width$legend$bottom$args <- x$legend$bottom$args
  ##   M.width$xlab.top <- list(K.ylab, just=1, col='transparent')
  ##   M.width$plot.args$panel.width <- NULL
  ##   print(M.width, position=c(px$M[1], 0, px$M[2], 1), more=FALSE)
  ## }
#  
#  
  print(L, position=c(px$L[1], 0, px$L[2], 1), more=TRUE)
  print(R, position=c(px$R[1], 0, px$R[2], 1), more=TRUE)
  print(M, position=c(px$M[1], 0, px$M[2], 1), more=FALSE)

  invisible(x)
  
}

emptyLegend <- function(args) { ## empty a legend and keep the space allocated to the legend.
  if (is.null(args)) args
  args$text[] <- " "
  args$rect <- list(border=0)
  args$col <- 0
  args$rectangles <- FALSE
  if (!is.null(args$title))
    args$title <- " "
  args
}



if (FALSE) {
c49 <- USAge.table[75:1, 2:1, "1949"]/1000000

PL <- likert(c49,
             main="Population of United States 1979 (ages 0-74)",
             xlab="Count in Millions",
             ylab="Age",
             scales=list(
               y=list(
                 limits=c(0,77),
                 at=seq(1,76,5),
                 labels=seq(0,75,5),
                 tck=.5))
             )
PL
print.pyramidLikert(PL)
as.pyramidLikert(PL)
PL0 <- as.pyramidLikert(PL)
PL0
print(as.pyramidLikert(PL), panel.width=.45)
plot(PL0)

PL.no.ylab <- likert(c49,
                     main="Population of United States 1979 (ages 0-74)",
                     xlab="Count in Millions",
                     ylab="",
                     scales=list(
                       y=list(
                         limits=c(0,77),
                         at=seq(1,76,5),
                         labels=seq(0,75,5),
                       tck=.5))
                     )
PL.no.ylab
as.pyramidLikert(PL.no.ylab)

## alternate spacing
as.pyramidLikert(PL, panel.width=.4) ## too much space in middle
as.pyramidLikert(update(PL, xlim=c(-2.05, 2.05)), panel.width=.4) ## too much space in middle
as.pyramidLikert(update(PL, xlim=c(-2.05, 2.05)), panel.width=.6) ## silly with overlap

}

## source("c:/HOME/rmh/HH-R.package/HH/R/pyramidLikert.R")
