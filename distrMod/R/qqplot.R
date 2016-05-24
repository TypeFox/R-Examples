################################################################
# QQ - Plot functions in package distrMod
################################################################

### to be written into the respective MASKING files....


## helper into distrMod
.labelprep <- function(x,y,lab.pts,col.lbl,cex.lbl,which.lbs,which.Order,order.traf){
      n <- length(x)
      rx <- rank(x)
      xys <- cbind(x,y[rx])
      if(is.null(which.lbs)) which.lbs <- 1:n
      oN0 <- order(x,decreasing=TRUE)
      if(!is.null(order.traf)){
          oN0 <- order(order.traf(x),decreasing=TRUE)
      }
      oN0b <- oN0 %in% which.lbs
      oN0 <- oN0[oN0b]
      oN <- oN0
      if(!is.null(which.Order))
          oN <- oN0[which.Order]
      x0 <- xys[oN,1]
      y0 <- xys[oN,2]

      col.lbl <- col.lbl[rx]
      lab.pts <- lab.pts[rx]
      cex.lbl <- cex.lbl[rx]
      return(list(x0=x0,y0=y0,lab=lab.pts[oN],col=col.lbl[oN],cex=cex.lbl[oN]))
}


### helper functions from distr

.confqq          <- distr:::.confqq
.isReplicated    <- distr:::.isReplicated
.makeLenAndOrder <- distr:::.makeLenAndOrder
.inGaps          <- distr:::.inGaps
.deleteItemsMCL  <- distr:::.deleteItemsMCL
.NotInSupport    <- distr:::.NotInSupport

setMethod("qqplot", signature(x = "ANY",
                              y = "UnivariateDistribution"),
    function(x,    ### observations
             y,    ### distribution
             n = length(x), ### number of points to be plotted
             withIdLine = TRUE, ### shall line y=x be plotted in
             withConf = TRUE,   ### shall confidence lines be plotted
             withConf.pw  = withConf,   ### shall pointwise confidence lines be plotted
             withConf.sim = withConf,   ### shall simultaneous confidence lines be plotted
             plot.it = TRUE,    ### shall be plotted at all (inherited from stats::qqplot)
             xlab = deparse(substitute(x)), ## x-label
             ylab = deparse(substitute(y)), ## y-label
             ...,                 ## further parameters
             width = 10,          ## width (in inches) of the graphics device opened
             height = 5.5,        ## height (in inches) of the graphics device opened}
             withSweave = getdistrOption("withSweave"), ## logical: if \code{TRUE}
             ##               (for working with \command{Sweave}) no extra device is opened and height/width are not set
             mfColRow = TRUE,     ## shall we use panel partition mfrow=c(1,1)?
             n.CI = n,            ## number of points to be used for CI
             withLab = FALSE,     ## shall observation labels be plotted in
             lab.pts = NULL,      ## observation labels to be used
             which.lbs = NULL,    ## which observations shall be labelled
             which.Order = NULL,  ## which of the ordered (remaining) observations shall be labelled
             order.traf = NULL,   ## an optional trafo; by which the observations are ordered (as order(trafo(obs))
             col.IdL = "red",     ## color for the identity line
             lty.IdL = 2,         ## line type for the identity line
             lwd.IdL = 2,         ## line width for the identity line
             alpha.CI = .95,      ## confidence level
             exact.pCI = (n<100), ## shall pointwise CIs be determined with exact Binomial distribution?
             exact.sCI = (n<100), ## shall simultaneous CIs be determined with exact kolmogorov distribution?
             nosym.pCI = FALSE,   ## shall we use (shortest) asymmetric CIs?
             col.pCI = "orange",  ## color for the pointwise CI
             lty.pCI = 3,         ## line type for the pointwise CI
             lwd.pCI = 2,         ## line width for the pointwise CI
             pch.pCI = par("pch"),## symbol for points (for discrete mass points) in pointwise CI
             cex.pCI = par("cex"),## magnification factor for points (for discrete mass points) in pointwise CI
             col.sCI = "tomato2", ## color for the simultaneous CI
             lty.sCI = 4,         ## line type for the simultaneous CI
             lwd.sCI = 2,         ## line width for the simultaneous CI
             pch.sCI = par("pch"),## symbol for points (for discrete mass points) in simultaneous CI
             cex.sCI = par("cex"),## magnification factor for points (for discrete mass points) in simultaneous CI
             cex.pch = par("cex"),## magnification factor for the plotted symbols
             col.pch = par("col"),## color for the plotted symbols
             cex.lbl = par("cex"),## magnification factor for the plotted observation labels
             col.lbl = par("col"),## color for the plotted observation labels
             adj.lbl = NULL,      ## adj parameter for the plotted observation labels
             alpha.trsp = NA,     ## alpha transparency to be added afterwards
             jit.fac = 0,         ## jittering factor used for discrete distributions
             check.NotInSupport = TRUE, ## shall we check if all x lie in support(y)
             col.NotInSupport = "red", ## if preceding check TRUE color of x if not in support(y)
             with.legend = TRUE,  ## shall a legend be plotted
             legend.bg = "white", ## background for the legend
             legend.pos = "topleft", ## position for the legend
             legend.cex = 0.8,     ## magnification factor for the legend
             legend.pref = "",     ## prefix for legend  text
             legend.postf = "",    ## postfix for legend text
             legend.alpha = alpha.CI ## nominal level of CI
    ){ ## return value as in stats::qqplot

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]
    mcl$withSweave <- NULL
    mcl$mfColRow <- NULL
    force(x)


    xj <- x
    if(any(.isReplicated(x)))
       xj[.isReplicated(x)] <- jitter(x[.isReplicated(x)], factor=jit.fac)

    ord.x <- order(xj)

    pp <- ppoints(n)
    yc <- q(y)(pp)

    yc.o <- yc

    if("support" %in% names(getSlots(class(y))))
       yc <- sort(jitter(yc, factor=jit.fac))

    alp.v <- .makeLenAndOrder(alpha.trsp,ord.x)
    alp.t <- function(x,a1) if(is.na(x)) x else addAlphTrsp2col(x,a1)
    alp.f <- if(length(alpha.trsp)==1L && is.na(alpha.trsp))
             function(x,a) x else function(x,a) mapply(x,alp.t,a1=a)
    cex.pch <- .makeLenAndOrder(cex.pch,ord.x)
    cex.lbl <- .makeLenAndOrder(cex.lbl,ord.x)
    col.pch <- alp.f(.makeLenAndOrder(col.pch,ord.x),alp.v)
    col.lbl <- alp.f(.makeLenAndOrder(col.lbl,ord.x),alp.v)

    if(withLab){
      if(is.null(lab.pts)) lab.pts <- paste(ord.x)
      else lab.pts <- .makeLenAndOrder(lab.pts,ord.x)
    }

    if(check.NotInSupport){
       xo <- x[ord.x]
       nInSupp <- which(xo < q(y)(0))

       nInSupp <- unique(sort(c(nInSupp,which( xo > q(y)(1)))))
       if("support" %in% names(getSlots(class(y))))
          nInSupp <- unique(sort(c(nInSupp,which( ! xo %in% support(y)))))
       if("gaps" %in% names(getSlots(class(y))))
          nInSupp <- unique(sort(c(nInSupp,which( .inGaps(xo,gaps(y))))))
       if(length(nInSupp)){
          col.pch[nInSupp] <- col.NotInSupport
          if(withLab)
#             col.lbl[ord.x[nInSupp]] <- col.NotInSupport
             col.lbl[nInSupp] <- col.NotInSupport
       }
    }


    if(n!=length(x)) withLab <- FALSE

    mcl$x <- xj
    mcl$y <- yc
    mcl <- .deleteItemsMCL(mcl)
    mcl$cex <- cex.pch
    mcl$col <- col.pch

    if (!withSweave){
           devNew(width = width, height = height)
    }
    opar <- par("mfrow", no.readonly = TRUE)
    if(mfColRow) on.exit(do.call(par, list(mfrow=opar, no.readonly = TRUE)))

    if(mfColRow) opar1 <- par(mfrow = c(1,1), no.readonly = TRUE)

    ret <- do.call(stats::qqplot, args=mcl)

    if(withLab&& plot.it){
       lbprep <- .labelprep(xj,yc,lab.pts,
                            col.lbl,cex.lbl,which.lbs,which.Order,order.traf)
       text(x = lbprep$x0, y = lbprep$y0, labels = lbprep$lab,
            cex = lbprep$cex, col = lbprep$col, adj = adj.lbl)
    }

    if(withIdLine&& plot.it){
       abline(0,1,col=col.IdL,lty=lty.IdL,lwd=lwd.IdL)
       if(#is(y,"AbscontDistribution")&&
       withConf){
          xy <- unique(sort(c(x,yc.o)))
          xy <- xy[!.NotInSupport(xy,y)]
          lxy <- length(xy)
          if(is(y,"DiscreteDistribution")){
             n0 <- min(n.CI, length(support(y)))
             n1 <- max(n0-lxy,0)
             if (n1 >0 ){
                 notyetInXY <- setdiff(support(y), xy)
                 xy0 <- sample(notyetInXY, n1)
                 xy <- sort(unique(c(xy,xy0)))
             }
          }else{
             if(lxy < n.CI){
                n1 <- (n.CI-lxy)%/%3
                xy0 <- seq(min(xy),max(xy),length=n1)
                xy1 <- r(y)(n.CI-lxy-n1)
                xy <- sort(unique(c(xy,xy0,xy1)))
             }
          }

          .confqq(xy, y, withConf.pw, withConf.sim, alpha.CI,
                      col.pCI, lty.pCI, lwd.pCI, pch.pCI, cex.pCI,
                      col.sCI, lty.sCI, lwd.sCI, pch.sCI, cex.sCI,
                  n, exact.sCI = exact.sCI, exact.pCI = exact.pCI,
                  nosym.pCI = nosym.pCI, with.legend = with.legend,
                  legend.bg = legend.bg, legend.pos = legend.pos,
                  legend.cex = legend.cex, legend.pref = legend.pref,
                  legend.postf = legend.postf, legend.alpha = legend.alpha)
       }
    }
    return(ret)
    })

## into distrMod
setMethod("qqplot", signature(x = "ANY",
                              y = "ProbFamily"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...){

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]

    mcl$y <- yD <- y@distribution
    if(!is(yD,"UnivariateDistribution"))
       stop("Not yet implemented.")

    return(do.call(getMethod("qqplot", signature(x="ANY", y="UnivariateDistribution")),
            args=mcl))
    })

