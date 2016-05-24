setMethod("plot", signature(x = "IC", y = "missing"),
    function(x, ...,withSweave = getdistrOption("withSweave"),
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.legend = FALSE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             withMBR = FALSE, MBRB = NA, MBR.fac = 2, col.MBR = par("col"),
             lty.MBR = "dashed", lwd.MBR = 0.8,
             scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL){

        xc <- match.call(call = sys.call(sys.parent(1)))$x
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
        dotsLeg <- dotsT <- dotsL <- .makedotsLowLevel(dots)

        if(!is.logical(inner)){
          if(!is.list(inner))
              inner <- as.list(inner)
            #stop("Argument 'inner' must either be 'logical' or a 'list'")
           inner <- .fillList(inner,4)
           innerD <- inner[1:3]
           innerL <- inner[4] 
        }else{innerD <- innerL <- inner}


        L2Fam <- eval(x@CallL2Fam)
        if(missing(scaleX.fct)){
           scaleX.fct <- p(L2Fam)
           scaleX.inv <- q(L2Fam)
        }

        trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO)
        
        to.draw <- 1:dims
        dimnms  <- c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg)) 
                 to.draw <- pmatch(to.draw.arg, dimnms)
            else if(is.numeric(to.draw.arg)) 
                 to.draw <- to.draw.arg
        }
        dims0 <- length(to.draw)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)

        if(!is.null(x.ticks)) dots$xaxt <- "n"
        if(!is.null(y.ticks)){
           y.ticks <- .fillList(list(y.ticks), dims0)
           dots$yaxt <- "n"
        }

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        MBRB <- MBRB * MBR.fac

        e1 <- L2Fam@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        if(is(e1, "UnivariateDistribution")){
           xlim <- eval(dots$xlim)
           if(!is.null(xlim)){ 
               xm <- min(xlim)
               xM <- max(xlim)
               if(!length(xlim) %in% c(2,2*dims0))
                  stop("Wrong length of Argument xlim");
               xlim <- matrix(xlim, 2,dims0)
            }
            if(is(e1, "AbscontDistribution")){
                lower0 <- getLow(e1, eps = getdistrOption("TruncQuantile")*2)
                upper0 <- getUp(e1, eps = getdistrOption("TruncQuantile")*2)
                me <- median(e1); s <- IQR(e1)
                lower1 <- me - 6 * s
                upper1 <- me + 6 * s
                lower <- max(lower0, lower1)
                upper <- min(upper0, upper1)
                if(!is.null(xlim)){ 
                  lower <- min(lower,xm)
                  upper <- max(upper,xM)
                }
                h <- upper - lower
                x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
                plty <- "l"
                lty <- "solid"
            }else{
                if(is(e1, "DiscreteDistribution")) x.vec <- support(e1)
                else{
                   x.vec <- r(e1)(1000)
                   x.vec <- sort(unique(x.vec))
                }
                plty <- "p"
                lty <- "dotted"
                if(!is.null(dots$xlim)) x.vec <- x.vec[(x.vec>=xm) & (x.vec<=xM)]

            }
         }
         ylim <- eval(dots$ylim)
         if(!is.null(ylim)){ 
               if(!length(ylim) %in% c(2,2*dims0)) 
                  stop("Wrong length of Argument ylim"); 
               ylim <- matrix(ylim, 2,dims0)
         }

        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        xlab <- dots$xlab; if(is.null(xlab)) xlab <- "x"
        ylab <- dots$ylab; if(is.null(ylab)) ylab <- "(partial) IC"
        dots$xlab <- dots$ylab <- NULL

        IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")

        mainL <- FALSE
        subL <- FALSE
        lineT <- NA

     .mpresubs <- function(inx)
                    .presubs(inx, c("%C", "%D", "%A"),
                          c(as.character(class(x)[1]),
                            as.character(date()),
                            as.character(deparse(xc))))

     if (hasArg(main)){
         mainL <- TRUE
         if (is.logical(main)){
             if (!main) mainL <-  FALSE
             else
                  main <- gettextf("Plot for IC %%A") ###
                          ### double  %% as % is special for gettextf
             }
         main <- .mpresubs(main)
         if (mainL) {
             if(missing(tmar))
                tmar <- 5
             if(missing(cex.inner))
                cex.inner <- .65
             lineT <- 0.6
             }
     }
     if (hasArg(sub)){
         subL <- TRUE
         if (is.logical(sub)){
             if (!sub) subL <-  FALSE
             else       sub <- gettextf("generated %%D")
                          ### double  %% as % is special for gettextf
         }
         sub <- .mpresubs(sub)
         if (subL)
             if (missing(bmar)) bmar <- 6
     }

     if(is.logical(innerL)){
        tnm  <- c(rownames(trafO))
        tnms <- if(is.null(tnm)) paste(1:dims) else 
                                 paste("'", tnm, "'", sep = "") 
        mnm <- names(L2Fam@param@main)
        mnms <- if(is.null(mnm)) NULL else paste("'", mnm, "' = ", sep = "") 
        mss  <- paste(mnms, round(L2Fam@param@main, 3), collapse=", ",sep="")
        innerT <- paste(gettextf("Component "),  tnms, 
                        gettextf("\nof"), #gettextf(" of L_2 derivative\nof"),
                        name(x)[1],
                        gettextf("\nwith main parameter ("), mss,")")
        if(!is.null(L2Fam@param@nuisance)){
            nnm <- names(L2Fam@param@nuisance)
            nnms <- if(is.null(nnm)) NULL else paste("'", nnm, "' = ", sep = "") 
            innerT <- paste(innerT,
                        gettextf("\nand nuisance parameter ("),
                        paste(nnms,round(L2Fam@param@nuisance, 3), collapse = ", "),
                        ")",
                        sep=""  )
        }
        if(!is.null(L2Fam@param@fixed)){
            fnm <- names(L2Fam@param@fixed)
            fnms <- if(is.null(fnm)) NULL else paste("'", fnm, "' = ", sep = "") 
            innerT <- paste(innerT,
                        gettextf("\nand fixed known parameter ("),
                        paste(fnms, round(L2Fam@param@fixed, 3), collapse = ", "),
                        ")",
                        sep=""  )
        }
     }else{
        innerT <- lapply(inner, .mpresubs)
        innerT <- .fillList(innerT,dims)
        if(dims0<dims){
           innerT0 <- innerT
           for(i in 1:dims0) innerT[to.draw[i]] <- innerT0[i]          
        }
     }

        if(with.legend){
          fac.leg <- if(dims0>1) 3/4 else .75/.8
          if(missing(legend.location)){
             legend.location <- .fillList(list("bottomright"), dims0)
          }else{
             legend.location <- as.list(legend.location)
             legend.location <- .fillList(legend.location, dims0)
          }
          if(is.null(legend)){
             legend <- vector("list",dims0)
             legend <- .fillList(as.list(xc),dims0)
          }
        }


        w0 <- getOption("warn")
        options(warn = -1)
        on.exit(options(warn = w0))
        opar <- par(no.readonly = TRUE)
#        opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
        on.exit(par(opar))
        if (!withSweave)
             devNew()
        
        parArgs <- NULL
        if(mfColRow)
           parArgs <- list(mfrow = c(nrows, ncols))

        omar <- par("mar")
        parArgs <- c(parArgs,list(mar = c(bmar,omar[2],tmar,omar[4])))

        do.call(par,args=parArgs)


        dotsT["pch"] <- dotsT["cex"] <- NULL
        dotsT["col"] <- dotsT["lwd"] <- NULL
        dotsL["cex"] <- dotsLeg["bg"] <- dotsLeg["cex"] <- NULL
        dots$ylim <- NULL

        for(i in 1:dims0){
            indi <- to.draw[i]
            if(!is.null(ylim)) dots$ylim <- ylim[,i]       
            fct <- function(x) sapply(x, IC1@Map[[indi]])
            print(xlim[,i])
            resc <-.rescalefct(x.vec, fct, scaleX, scaleX.fct,
                              scaleX.inv, scaleY, scaleY.fct, xlim[,i],
                              ylim[,i], dots)
            dots <- resc$dots
            dots$xlim <- xlim[,i]
            dots$ylim <- ylim[,i]
            x.vec1 <- resc$X
            y.vec1 <- resc$Y
            do.call(plot, args=c(list(x=x.vec1, y=y.vec1, type = plty, lty = lty,
                                      xlab = xlab, ylab = ylab), dots))
            .plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct, scaleY.inv,
                              xlim[,i], ylim[,i], x.vec1, ypts = 400, n = scaleN,
                              x.ticks = x.ticks, y.ticks = y.ticks[[i]])
            if(withMBR){
                MBR.i <- MBRB[i,]
                if(scaleY) MBR.i <- scaleY.fct(MBR.i)
                abline(h=MBR.i, col=col.MBR, lty=lty.MBR, lwd = lwd.MBR)
            }
            if(is(e1, "DiscreteDistribution")){
                x.vec1D <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
                rescD <-.rescalefct(x.vec1D, fct, scaleX, scaleX.fct,
                                scaleX.inv, scaleY, scaleY.fct, xlim[,i],
                                ylim[,i], dots)
                x.vecD <- rescD$X
                y.vecD <- rescD$Y

                dotsL$lty <- NULL
                do.call(lines,args=c(list(x.vecD, y.vecD,
                                          lty = "dotted"), dotsL))
            }
            do.call(title,args=c(list(main = innerT[indi]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner))
            if(with.legend)
               legend(.legendCoord(legend.location[[i]], scaleX, scaleX.fct,
                        scaleY, scaleY.fct), bg = legend.bg,
                      legend = legend[[i]], dotsLeg, cex = legend.cex*fac.leg)

        }
        cex.main <- if(!hasArg(cex.main)) par("cex.main") else dots$"cex.main"
        col.main <- if(!hasArg(col.main)) par("col.main") else dots$"col.main"
        if (mainL)
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)

        cex.sub <- if(!hasArg(cex.sub)) par("cex.sub") else dots$"cex.sub"
        col.sub <- if(!hasArg(col.sub)) par("col.sub") else dots$"col.sub"
        if (subL)
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)

        invisible()
    })


setMethod("plot", signature(x = "IC",y = "numeric"),
          function(x, y, ..., cex.pts = 1, col.pts = par("col"),
          pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
          lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
          which.lbs = NULL, which.Order  = NULL, return.Order = FALSE){

    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."

    n <- if(!is.null(dim(y))) nrow(y) else length(y)
    pch.pts <- rep(pch.pts, length.out=n)
    lab.pts <- if(is.null(lab.pts)) paste(1:n) else rep(lab.pts,n)


    L2Fam <- eval(x@CallL2Fam)
    trafO <- trafo(L2Fam@param)
    dims <- nrow(trafO)
    dimm <- length(L2Fam@param)
    QF <- diag(dims)

    if(is(x,"ContIC") & dims>1 )
      {if (is(normtype(x),"QFNorm")) QF <- QuadForm(normtype(x))}

    IC1 <- as(diag(dims) %*% x@Curve, "EuclRandVariable")
    absInfo <- t(IC1) %*% QF %*% IC1
    ICMap <- IC1@Map

    sel <- .SelectOrderData(y, function(x)sapply(x, absInfo@Map[[1]]),
                            which.lbs, which.Order)
    i.d <- sel$ind
    i0.d <- sel$ind1
    n <- length(i.d)

    dots.without <- dots
    dots.without$col <- dots.without$cex <- dots.without$pch <- NULL


    pL <- expression({})
    if(!is.null(dots$panel.last))
        pL <- dots$panel.last
    dots$panel.last <- NULL

    pL <- substitute({
        y1 <- y0s
        ICy <- sapply(y0s,ICMap0[[indi]])
        print(xlim[,i])
        resc.dat <-.rescalefct(y0s, function(x) sapply(x,ICMap0[[indi]]),
                              scaleX, scaleX.fct, scaleX.inv,
                              scaleY, scaleY.fct, xlim[,i], ylim[,i],
                              dwo0)
        y1 <- resc.dat$X
        ICy <- resc.dat$Y

        if(is(e1, "DiscreteDistribution"))
           ICy <- jitter(ICy, factor = jitter.fac0)

        col.pts <- if(!is.na(al0)) sapply(col0, addAlphTrsp2col,alpha=al0) else col0

        do.call(points, args=c(list(y1, ICy, cex = log(absy0+1)*3*cex0,
                        col = col.pts, pch = pch0), dwo0))
        if(with.lab0){
           text(x = y0s, y = ICy, labels = lab.pts0,
                cex = log(absy0+1)*1.5*cex0, col = col0)
        }
        pL0
        }, list(pL0 = pL, ICMap0 = ICMap, y0s = sel$data, absy0 = sel$y,
                dwo0 = dots.without, cex0 = cex.pts, pch0 = pch.pts[i.d],
                col0 = col.pts, with.lab0 = with.lab, lab.pts0 = lab.pts[i.d],
                al0 = alpha.trsp, jitter.fac0 = jitter.fac
                ))

  do.call("plot", args = c(list(x = x, panel.last = pL), dots))
  if(return.Order) return(i0.d)
  invisible()
})

