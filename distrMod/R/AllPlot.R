### from Matthias' thesis / ROptEst
setMethod("plot", signature(x = "ParamFamily", y = "missing"),
    function(x, ...){ 
        e1 <- x@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        plot(e1) 
    })
setMethod("plot", signature(x = "L2ParamFamily", y = "missing"),
    function(x, withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], ...,
             mfColRow = TRUE, to.draw.arg = NULL){

        xc <- match.call(call = sys.call(sys.parent(1)))$x
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
        
        dots$to.draw.arg <- NULL
        trafO <- trafo(x@param)
#        dims <- nrow(trafO)
        dimm <- dims <- length(x@param)
        
        to.draw <- 1:(3+dims)
        dimnms  <- names(main(x@param)) #c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        names(to.draw) <- c("d","p","q", dimnms)
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg)) 
                 to.draw <- pmatch(to.draw.arg, names(to.draw))
            else if(is.numeric(to.draw.arg)) 
                 to.draw <- to.draw.arg
        }
        l2dpl <- to.draw[to.draw > 3]
        dims0 <- length(l2dpl)
        nrows <- trunc(sqrt(dims0))
        ncols <- ceiling(dims0/nrows)

        
        if(!is.logical(inner)){
          if(!is.list(inner))
              inner <-  as.list(inner)
            #stop("Argument 'inner' must either be 'logical' or a 'list'")
          innerLog <- TRUE  
          iL <- length(to.draw[to.draw <= 3])+length(l2dpl)
          iLD <- (1:iL)[to.draw <= 3]
          iLL <- (1:iL)[to.draw > 3]
          inner <- distr:::.fillList(inner,iL)          
          innerD <- if(length(iLD)) inner[iLD] else NULL
          innerL <- if(length(iLL)) inner[iLL] else NULL
        }else{innerLog <- innerD <- innerL <- inner}
        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL

        e1 <- x@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")


        if(is(e1, "UnivariateDistribution")){
           xlim <- eval(dots$xlim)
           if(!is.null(xlim)){ 
               xm <- min(xlim)
               xM <- max(xlim)
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
                if(!is.null(xlim)) x.vec <- x.vec[(x.vec>=xm) & (x.vec<=xM)]
            }
        }
        dxg <- d(e1)(x.vec)
        pxg <- p(e1)(x.vec)
        ylim <- eval(dots$ylim)
        if(!is.null(ylim)){ 
               d.0 <- 1 %in% to.draw
               d.1 <- 2 %in% to.draw | 3 %in% to.draw
               if(! length(ylim) %in% c(2,2*(d.0+d.1+dims0))) 
                  stop("Wrong length of Argument ylim"); 
               ylim <- matrix(ylim, 2,d.0+d.1+dims0)
               iy <- if(d.0+d.1==2) 1:2 else 1
               dots$ylim <- ylim[,iy]
        }

        
        L2deriv <- as(diag(dimm) %*% x@L2deriv, "EuclRandVariable")

        mainL <- FALSE
        subL <- FALSE
        lineT <- NA

     .mpresubs <- function(inx)
                    distr:::.presubs(inx, c("%C", "%D", "%A"),
                          c(as.character(class(x)[1]),
                            as.character(date()),
                            as.character(deparse(xc))))

     if (hasArg(main)){
         mainL <- TRUE
         if (is.logical(main)){
             if (!main) mainL <-  FALSE
             else
                  main <- gettextf("Distribution Plot for model %%A") ###
                          ### double  %% as % is special for gettextf
             }
         main <- .mpresubs(main)
         if (mainL) {
             if(missing(tmar))
                tmar <- 5
             lineT <- 0.6
             }
     }
     if(missing(cex.inner)){
        cex.inner <- .65
        cex.innerD <- 1
     }else{
        cex.inner <- rep(cex.inner, length.out=2)
        cex.innerD <- cex.inner[1]
        cex.inner <- cex.inner[2]             
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
#        tnm  <- c(rownames(trafO))
        tnms <- 1:dims #if(is.null(tnm)) paste(1:dims) else paste("'", tnm, "'", sep = "") 
        mnm <- names(x@param@main)
        mnms <- if(is.null(mnm)) NULL else paste("'", mnm, "' = ", sep = "") 
        mss  <- paste(mnms, round(x@param@main, 3), collapse=", ",sep="")
        innerT <- paste(gettextf("Component "),  tnms, 
                        gettextf(" of L_2 derivative\nof"),
                        name(x)[1],
                        gettextf("\nwith main parameter ("), mss,")")
        if(!is.null(x@param@nuisance)){
            nnm <- names(x@param@nuisance)
            nnms <- if(is.null(nnm)) NULL else paste("'", nnm, "' = ", sep = "") 
            innerT <- paste(innerT,
                        gettextf("\nand nuisance parameter ("),
                        paste(nnms,round(x@param@nuisance, 3), collapse = ", "),
                        ")",
                        sep=""  )
        }
        if(!is.null(x@param@fixed)){
            fnm <- names(x@param@fixed)
            fnms <- if(is.null(fnm)) NULL else paste("'", fnm, "' = ", sep = "") 
            innerT <- paste(innerT,
                        gettextf("\nand fixed known parameter ("),
                        paste(fnms, round(x@param@fixed, 3), collapse = ", "),
                        ")",
                        sep=""  )
        }
        innerT <- if(length(l2dpl)) innerT[l2dpl-3] else NULL
     }else{
        innerT <- lapply(innerL, .mpresubs)
        innerD <- lapply(innerD, .mpresubs)
     }


        dotsT <- dots
        dotsT["main"] <- NULL
        dotsT["cex.main"] <- NULL
        dotsT["col.main"] <- NULL
        dotsT["line"] <- NULL

        distrpl <- (1:3) %in% to.draw
        todrw <- as.numeric((1:3)[distrpl])
        if(any(distrpl)){
           lis0 <- c(list(e1, withSweave = withSweave, 
                          main = main, inner = innerD, sub = sub, 
                          col.inner = col.inner, cex.inner = cex.innerD),
                     dots, mfColRow = mfColRow)
           lis0$to.draw.arg  <- todrw 
           do.call(plot, args = lis0)            
        }
        o.warn <- options("warn")
        options(warn = -1)
        on.exit(options(warn=o.warn))
        opar <- par(no.readonly = TRUE)
   #     opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
        on.exit(par(opar, no.readonly = TRUE))
        
        if (!withSweave)
             devNew()
        
        parArgs <- NULL
        if(mfColRow)
           parArgs <- list(mfrow = c(nrows, ncols))

        omar <- par("mar", no.readonly = TRUE)
        parArgs <- c(parArgs,list(mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE))
       
        dots$ylim <- NULL
        do.call(par,args=parArgs)
        for(i in 1:dims0){
            indi <- l2dpl[i]-3
            if(!is.null(ylim)) dots$ylim <- ylim[,d.0+d.1+i]       
            do.call(plot, args=c(list(x=x.vec, y=sapply(x.vec, L2deriv@Map[[indi]]),
                                 type = plty, lty = lty,
                                 xlab = "x",
                                 ylab = expression(paste(L[2], " derivative"))),
                                 dots))
            if(is(e1, "DiscreteDistribution")){
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
                do.call(lines, args=c(list(x.vec1, sapply(x.vec1, L2deriv@Map[[indi]]),
                              lty = "dotted"),dots))
            }
            if(innerLog)
               do.call(title, args = c(list(main = innerT[i]), dotsT, 
                       line = lineT, cex.main = cex.inner, 
                       col.main = col.inner))
        }

        if(!hasArg(cex.main)) cex.main <- par("cex.main") else cex.main <- dots$"cex.main"
        if(!hasArg(col.main)) col.main <- par("col.main") else col.main <- dots$"col.main"
        if (mainL)
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)

        if(!hasArg(cex.sub)) cex.sub <- par("cex.sub") else cex.sub <- dots$"cex.sub"
        if(!hasArg(col.sub)) col.sub <- par("col.sub") else col.sub <- dots$"col.sub"
        if (subL)
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)

     invisible()
    })
