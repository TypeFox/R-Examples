################################################################
# QQ - Plot functions in package distr
################################################################

setMethod("qqplot", signature(x = "UnivariateDistribution",
                              y = "UnivariateDistribution"), function(x, y,
                              n = 30, withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...,
    width = 10, height = 5.5, withSweave = getdistrOption("withSweave"),
    mfColRow=TRUE, n.CI = n,
    col.IdL = "red", lty.IdL = 2, lwd.IdL = 2,
    alpha.CI = .95, exact.pCI = (n<100), exact.sCI = (n<100), nosym.pCI = FALSE,
    col.pCI = "orange", lty.pCI = 3, lwd.pCI = 2, pch.pCI = par("pch"), cex.pCI = par("cex"),
    col.sCI = "tomato2", lty.sCI = 4, lwd.sCI = 2, pch.sCI = par("pch"), cex.sCI = par("cex"),
    cex.pch = par("cex"), col.pch = par("col"),
    jit.fac = 0, check.NotInSupport = TRUE,
    col.NotInSupport = "red", with.legend = TRUE, legend.bg = "white",
    legend.pos = "topleft", legend.cex = 0.8, legend.pref = "", 
    legend.postf = "", legend.alpha = alpha.CI, debug = FALSE, withSubst = TRUE){

    mc <- match.call(call = sys.call(sys.parent(1)))
    xcc <- as.character(deparse(mc$x))
    ycc <- as.character(deparse(mc$y))
    if(missing(xlab)) mc$xlab <- xcc
    if(missing(ylab)) mc$ylab <- ycc

    mcl <- as.list(mc)[-1]
    mcl$withSweave <- NULL
    mcl$mfColRow <- NULL
    mcl$debug <- NULL

   .mpresubs <- if(withSubst){
                 function(inx) 
                    .presubs(inx, c("%C", "%A", "%D" ),
                          c(as.character(class(x)[1]), 
                            as.character(date()), 
                            xcc))
                }else function(inx) inx

    force(x)

    pp <- ppoints(n)
    xc <- q(x)(pp)
    yc <- q(y)(pp)

    col.pch <- rep(col.pch,length.out=n)

    if(check.NotInSupport){
       xco <- sort(xc)
       nInSupp <- .NotInSupport(xc,y)
       if(length(nInSupp)){
          col.pch[nInSupp] <- col.NotInSupport
       }
    }


    oxc <- 1:length(xc)
    xc.o <- xc
    yc.o <- yc
    ord.x <- order(xc)

    if("support" %in% names(getSlots(class(x)))){
       xc <- jitter(xc, factor=jit.fac)
       oxc <- order(xc)
       xc <- xc[oxc]
       }

    if("support" %in% names(getSlots(class(y))))
       yc <- sort(jitter(yc, factor=jit.fac))

    mcl$x <- xc
    mcl$y <- yc

    mcl <- .deleteItemsMCL(mcl)
    mcl$cex <- .makeLenAndOrder(cex.pch,ord.x)
    mcl$col <- .makeLenAndOrder(col.pch,ord.x)

    mcl$xlab <- .mpresubs(mcl$xlab)
    mcl$ylab <- .mpresubs(mcl$ylab)

    if (!is.null(eval(mcl$main)))
        mcl$main <- .mpresubs(eval(mcl$main))
    if (!is.null(eval(mcl$sub)))
        mcl$sub <- .mpresubs(eval(mcl$sub))


    if (!withSweave){
           devNew(width = width, height = height)
    }
    opar <- par("mfrow", no.readonly = TRUE)
    if(mfColRow) on.exit(do.call(par, list(mfrow=opar)))
    if(mfColRow) opar1 <- par(mfrow = c(1,1), no.readonly = TRUE)

    ret <- do.call(stats::qqplot, args=mcl)
    qqb <- NULL
    if(withIdLine){
       if(plot.it)abline(0,1,col=col.IdL,lty=lty.IdL,lwd=lwd.IdL)
       if(#is(y,"AbscontDistribution") &&
       withConf){
          xy <- unique(sort(c(xc.o,yc.o)))
          xy <- xy[!.NotInSupport(xy,y)]
          lxy <- length(xy)
          if(is(y,"DiscreteDistribution")){
             n0 <- min(n.CI, length(support(y)))
             n1 <- max(n0-lxy,0)
             if (n1 >0 ){
                 notyetInXY <- setdiff(support(y), xy)
                 xy0 <- sample(notyetInXY, n1)
                 xy <- sort(c(xy,xy0))
             }
          }else{
             if(lxy < n.CI){
                n1 <- (n.CI-lxy)%/%3
                xy0 <- seq(min(xy),max(xy),length=n1)
                xy1 <- r(y)(n.CI-lxy-n1)
                xy <- sort(c(xy,xy0,xy1))
             }
          }
        qqplotInfo <- list(xy.0=xy, y.0=y, 
                         withConf.pw=withConf.pw, 
                         withConf.sim=withConf.sim, 
                         alpha.CI=alpha.CI ,
                         col.pCI = col.pCI , lty.pCI = lty.pCI , 
                         lwd.pCI = lwd.pCI , pch.pCI = pch.pCI, 
                         cex.pCI = cex.pCI , 
                         col.sCI = col.sCI , lty.sCI = lty.sCI , 
                         lwd.sCI = lwd.sCI , pch.sCI = pch.sCI, 
                         cex.sCI = cex.sCI , 
                         n = n , 
                         exact.sCI = exact.sCI, exact.pCI = exact.pCI,
                  nosym.pCI = nosym.pCI, with.legend = with.legend,
                  legend.bg = legend.bg, legend.pos = legend.pos,
                  legend.cex = legend.cex, legend.pref = legend.pref,
                  legend.postf = legend.postf, legend.alpha = legend.alpha, debug = debug,
                  args.stats.qqplot = mcl
                  )
       if(plot.it){
           qqb <- .confqq(xy, y, datax=TRUE, withConf.pw, withConf.sim, alpha.CI,
                      col.pCI, lty.pCI, lwd.pCI, pch.pCI, cex.pCI,
                      col.sCI, lty.sCI, lwd.sCI, pch.sCI, cex.sCI,
                  n, exact.sCI = exact.sCI, exact.pCI = exact.pCI,
                  nosym.pCI = nosym.pCI, with.legend = with.legend,
                  legend.bg = legend.bg, legend.pos = legend.pos,
                  legend.cex = legend.cex, legend.pref = legend.pref,
                  legend.postf = legend.postf, legend.alpha = legend.alpha,
                  debug = debug)
          }else{
           qqb <- qqbounds(sort(unique(xy)),y,alpha.CI,n,withConf.pw, withConf.sim,
                   exact.sCI,exact.pCI,nosym.pCI,debug)
          }
       }
    }
    qqplotInfo <- c(ret, qqplotInfo, qqb)
    class(qqplotInfo) <- c("qqplotInfo","DiagnInfo")
    return(invisible(qqplotInfo))
    })
    
