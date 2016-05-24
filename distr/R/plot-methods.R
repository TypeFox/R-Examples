# -------- AbscontDistribution ---------- #

setMethod("plot", signature(x = "AbscontDistribution", y = "missing"),
   function(x, width = 10, height = 5.5, withSweave = getdistrOption("withSweave"),
            xlim = NULL, ylim = NULL, ngrid = 1000, verticals = TRUE,
            do.points = TRUE, main = FALSE, inner = TRUE, sub = FALSE, 
            bmar = par("mar")[1], tmar = par("mar")[3], ..., 
            cex.main = par("cex.main"), cex.inner = 1.2, 
            cex.sub = par("cex.sub"), col.points = par("col"), 
            col.vert = par("col"), col.main = par("col.main"), 
            col.inner = par("col.main"), col.sub = par("col.sub"), 
            cex.points = 2.0, pch.u = 21, pch.a = 16, mfColRow = TRUE,
            to.draw.arg = NULL, withSubst = TRUE){

     xc <- match.call(call = sys.call(sys.parent(1)))$x
     ### manipulating the ... - argument
     dots <- match.call(call = sys.call(sys.parent(1)), 
                        expand.dots = FALSE)$"..."

      to.draw <- 1:3
      names(to.draw) <- c("d","p","q")
      if(! is.null(to.draw.arg)){
         if(is.character(to.draw.arg)) 
            to.draw <- pmatch(to.draw.arg, names(to.draw))
         else if(is.numeric(to.draw.arg)) 
              to.draw <- to.draw.arg
      }
      l.draw <- length(to.draw)

     
     pF <- expression({})
     if(!is.null(dots[["panel.first"]])){
         pF <- .panel.mingle(dots,"panel.first")
     }
     pF <- .fillList(pF, l.draw)
     pL <- expression({})
     if(!is.null(dots[["panel.last"]])){
          pL <- .panel.mingle(dots,"panel.last")
     }
     pL <- .fillList(pL, l.draw)

     dots$panel.first <- dots$panel.last <- NULL
     dots$col.hor <- NULL

     dots.for.points <- .makedotsPt(dots)

     dots.lowlevel <- .makedotsLowLevel(dots)
     dots.without.pch <- dots.lowlevel[! (names(dots.lowlevel) %in% c("col", "pch"))]
     if(!is(x,"AbscontDistribution"))
         x <- .ULC.cast(x)     
     ###
     if(!is.logical(inner))
         {if(!is.list(inner))
              inner <- as.list(inner)
            #stop("Argument 'inner' must either be 'logical' or a 'list'")
          inner <- .fillList(inner,l.draw)          
         }
     cex <- if (hasArg("cex")) dots$cex else 1

     if (hasArg("cex") && missing(cex.points)) 
         cex.points <- 2.0 * cex

     if (hasArg("pch") && missing(pch.u))
          pch.u <- dots$pch
     if (hasArg("pch") && missing(pch.a))
          pch.a <- dots$pch
     
     if (hasArg("col") && missing(col.points))
         col.points <- dots$col
     if (hasArg("col") && missing(col.vert))
         col.vert <- dots$col
     if (hasArg("col") && missing(col.main))
        col.main <- dots$col
     if (hasArg("col") && missing(col.inner))
        col.inner <- dots$col
     if (hasArg("col") && missing(col.sub))
        col.sub <- dots$col

     if (!withSweave){
           devNew(width = width, height = height)
           }
     omar <- par("mar", no.readonly = TRUE)
#     omar$cin <- omar$cra <- omar$csi <- omar$cxy <-  omar$din <- NULL
     if(mfColRow) on.exit(par(omar, no.readonly = TRUE))
     
     mainL <- FALSE
     subL <- FALSE
     lineT <- NA
     logpd <- logq <- ""

     slots <-  slotNames(param(x))
     slots <-  slots[slots != "name"]
     nrvalues <-  length(slots)
     if(nrvalues > 0){
           values <-  numeric(nrvalues)
       for(i in 1:nrvalues)
         values[i] <-  attributes(attributes(x)$param)[[slots[i]]]
       paramstring <-  paste(values, collapse = ", ")
       nparamstring <-  paste(slots, "=", values, collapse = ", ")    
       qparamstring <- paste("(",paramstring,")",sep="")
     }
     else paramstring <- qparamstring <- nparamstring <- ""

     .mpresubs <- if(withSubst){
             function(inx) 
                    .presubs(inx, c("%C", "%D", "%N", "%P", "%Q", "%A"),
                          c(as.character(class(x)[1]), 
                            as.character(date()), 
                            nparamstring, 
                            paramstring, 
                            qparamstring,
                            as.character(deparse(xc))))
            }else function(inx) inx
            
     xlab0 <- list("d"="x", "p"="q", "q"="p")
     iL <- 1:length(to.draw)
     .mp2 <- function(dlb = dots$xlab, lb0 = list("d"="x", "p"="q", "q"="p")){
              dlb0 <- eval(dlb)
              if (!is.null(dlb)){
              .mp <- if(is.list(dlb0)) function(x,i){
                                if(is.call(x)) x <- eval(x)
                                if(length(i)==0) return(NULL)
                                i <- min(i)
                                if(is.character(x[[i]])){
                                   return(as.character(eval(.mpresubs(x[[i]]))))
                                }else{
                                res <- .mpresubs(x[[i]])
                                if(length(res)==0) return(NULL)
                                if(is.call(res)) res <- res[-1]
                                return(res)}
                                }else function(x,i){
                                  if(length(x)==1) return(x[1])
                                  res <- x[i]
                                  if(length(res)==0) return(NULL)
                                  if(is.na(res)) return(NULL)
                                  return(res)}
              force(lb0)
              .mp3 <- .mp(dlb,iL[to.draw==1])
              if(1%in%to.draw & !is.null(.mp3)) lb0[["d"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==2])
              if(2%in%to.draw & !is.null(.mp3)) lb0[["p"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==3])
              if(3%in%to.draw & !is.null(.mp3)) lb0[["q"]] <- .mp3

             }
             return(lb0)}
     xlab0 <- .mp2()
     dots$xlab <- NULL
     ylab0 <- .mp2(dlb=dots$ylab, lb0=list("d"="d(x)", "p"="p(q)", "q"="q(p)"))
     dots$ylab <- NULL

     if (hasArg("main")){
         mainL <- TRUE
         if (is.logical(main)){
             if (!main) mainL <-  FALSE 
             else 
                  main <- gettextf("Distribution Plot for %%A") ### 
                          ### double  %% as % is special for gettextf
             }
         main <- .mpresubs(main)
         if (mainL) {
             if(missing(tmar))
                tmar <- 5
             if(missing(cex.inner)) 
                cex.inner <- .9 
             lineT <- 0.6
             }
     }
     if (hasArg("sub")){ 
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
     
     if(mfColRow)
         opar <- par(mfrow = c(1,l.draw), mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE)
     else
         opar <- par(mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE)
     
     if(is.logical(inner)){     
        inner.d <- if (inner) 
                   .mpresubs(gettextf("Density of %%C%%Q")) else ""
        inner.p <- if (inner) 
                   .mpresubs(gettextf("CDF of %%C%%Q")) else ""
        inner.q <- if (inner) 
                   .mpresubs(gettextf("Quantile function of %%C%%Q")) else ""
                          ### double  %% as % is special for gettextf
     }else{
        iL <- 1:length(to.draw)

        inner.d <- if(1%in%to.draw) .mpresubs(inner[[min(iL[to.draw==1])]]) else NULL
        inner.p <- if(2%in%to.draw) .mpresubs(inner[[min(iL[to.draw==2])]]) else NULL
        inner.q <- if(3%in%to.draw) .mpresubs(inner[[min(iL[to.draw==3])]]) else NULL
     }

     lower0 <- getLow(x, eps = getdistrOption("TruncQuantile")*2)
     upper0 <- getUp(x, eps = getdistrOption("TruncQuantile")*2)
     me <- q(x)(1/2); s <- q(x)(3/4)-q(x)(1/4)
     lower1 <- me - 6 * s
     upper1 <- me + 6 * s
     lower <- max(lower0, lower1)
     upper <- min(upper0, upper1)

     ## ngrid  nr of gridpoints
     ## exactq two p-values are considered equal if difference is
     ## is less than 10^-exactq in abs. value  

     dist <- upper - lower
     if(hasArg("xlim")) 
         {if(length(xlim)!=2) stop("Wrong length of Argument xlim");
             grid <- seq(xlim[1], xlim[2], length = ngrid)}
     else grid <- seq(from = lower - 0.1 * dist, to = upper + 0.1 * dist, 
                      length = ngrid)

     dxg <- d(x)(grid)
     pxg <- p(x)(grid)
     
     
     if(hasArg("ylim"))
         {if (3 %in% to.draw && any( c(1,2) %in% to.draw)){
                 if(! length(ylim) %in% c(2,4)) 
                     stop("Wrong length of Argument ylim")
           }else{
                 if(! length(ylim) == 2) 
                     stop("Wrong length of Argument ylim")
           }                  
           ylim <- matrix(ylim, 2,2)
           ylim1 <- ylim[,1]; ylim2 <- ylim[,2];
           }
     else {ylim1 <- c(0,max(dxg[dxg<50])); ylim2 <- c(-0.05,1.05)}

     if(hasArg("log"))
         {logpd <- dots$log
          logq <- gsub("u","y",gsub("y","x",gsub("x", "u", logpd)))
          if(length(grep("y",logpd))){ 
             ylim1 <- c(max(min(dxg[dxg>0]), ylim1[1]), 
                             ylim1[2])
             ylim2 <- c(max(min(pxg[pxg>0]), ylim2[1]), 
                             ylim2[2])
             }
          }

     plotCount <- 1
     o.warn <- getOption("warn"); options(warn = -1)
     if(1%in%to.draw){
         on.exit(options(warn=o.warn))
         dots.lowlevel$panel.first <- pF[[plotCount]]
         dots.lowlevel$panel.last  <- pL[[plotCount]]
         dots.lowlevel$xlim <- xlim
         do.call(plot, c(list(x = grid, dxg, type = "l",
             ylim = ylim1,  ylab = ylab0[["d"]], xlab = xlab0[["d"]], log = logpd),
             dots.lowlevel))
         dots.lowlevel$panel.first <- dots.lowlevel$panel.last <- NULL
         dots.lowlevel$xlim <- NULL
         plotCount <- plotCount + 1
         options(warn = o.warn)
     
         title(main = inner.d, line = lineT, cex.main = cex.inner,
               col.main = col.inner)
     
     
         options(warn = -1)
     }
     if(is.finite(q(x)(0))) {grid <- c(q(x)(0),grid); pxg <- c(0,pxg)}
     if(is.finite(q(x)(1))) {grid <- c(grid,q(x)(1)); pxg <- c(pxg,1)}

     if(2%in%to.draw){
        dots.lowlevel$panel.first <- pF[[plotCount]]
        dots.lowlevel$panel.last  <- pL[[plotCount]]
        dots.lowlevel$xlim <- xlim
        do.call(plot, c(list(x = grid, pxg, type = "l",
             ylim = ylim2, ylab = ylab0[["p"]], xlab = xlab0[["p"]], log = logpd),
             dots.lowlevel))
        dots.lowlevel$panel.first <- dots.lowlevel$panel.last <- NULL
        dots.lowlevel$xlim <- NULL
        plotCount <- plotCount + 1
        options(warn = o.warn)
      
        title(main = inner.p, line = lineT, cex.main = cex.inner,
              col.main = col.inner)
     }
     ### quantiles

     ### fix finite support bounds
     ixg  <-  grid>=max(q(x)(0),lower) & grid <= min(q(x)(1),upper)
     pxg  <-   pxg[ixg]
     grid <-  grid[ixg]

     ### fix constancy regions of p(x)
     if(isOldVersion(x)) x <- conv2NewVersion(x)
                 
     if(!is.null(gaps(x))){
        i.not.gap    <- !.isIn(grid,gaps(x))
        ndots <- nrow(gaps(x))
        pu1 <- p(x)(gaps(x)[,1])
        if (verticals){
             xu <- c(gaps(x)[,1],gaps(x)[,2], grid[i.not.gap])
             pu <- c(rep(pu1,2), pxg[i.not.gap])
        }else{
             xu <- c(gaps(x)[,1],rep(NA,ndots),gaps(x)[,2], grid[i.not.gap])
             pu <- c(rep(pu1,3), pxg[i.not.gap])
        }     
        #
        o <- order(pu)
        po <- pu[o]
        xo <- xu[o]
     }else{
        po <- pxg
        xo <- grid
     }
     
     if(3%in%to.draw){
        options(warn = -1)
        dots.lowlevel$panel.first <- pF[[plotCount]]
        dots.lowlevel$panel.last  <- pL[[plotCount]]
        do.call(plot, c(list(x = po, xo, type = "n",
             xlim = ylim2, ylim = xlim, ylab = ylab0[["q"]], xlab = xlab0[["q"]],
             log = logq), dots.lowlevel))
        dots.lowlevel$panel.first <- dots.lowlevel$panel.last <- NULL
        plotCount <- plotCount + 1
        options(warn = o.warn)
    
        
        title(main = inner.q, line = lineT, cex.main = cex.inner,
              col.main = col.inner)
        
        options(warn = -1)
            lines(po,xo, ...)
        if (verticals && !is.null(gaps(x))){
            pu <- rep(pu1,3)
            xu <- c(gaps(x)[,1],gaps(x)[,2],rep(NA,ndots))
            o <- order(pu)
            dots.without.pch0 <- dots.without.pch
            dots.without.pch0$col <- NULL
            do.call(lines, c(list(pu[o], xu[o], 
                    col = col.vert), dots.without.pch0))    
        }
        options(warn = o.warn)
     
        if(!is.null(gaps(x)) && do.points){
           do.call(points, c(list(x = pu1, y = gaps(x)[,1], pch = pch.a, 
                   cex = cex.points, col = col.points), dots.for.points) )
           do.call(points, c(list(x = pu1, y = gaps(x)[,2], pch = pch.u,
                   cex = cex.points, col = col.points), dots.for.points) )
        }
     }      
     if (mainL)
         mtext(text = main, side = 3, cex = cex.main, adj = .5, 
               outer = TRUE, padj = 1.4, col = col.main)                            
    
     if (subL)
         mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
               outer = TRUE, line = -1.6, col = col.sub)                            
   return(invisible())
   }
   )
# -------- DiscreteDistribution -------- #

setMethod("plot", signature(x = "DiscreteDistribution", y = "missing"),
    function(x,  width = 10, height = 5.5, withSweave = getdistrOption("withSweave"), 
             xlim = NULL, ylim = NULL, verticals = TRUE, do.points = TRUE, 
             main = FALSE, inner = TRUE, sub = FALSE,
             bmar = par("mar")[1], tmar = par("mar")[3], ..., 
             cex.main = par("cex.main"), cex.inner = 1.2, 
             cex.sub = par("cex.sub"), col.points = par("col"), 
             col.hor = par("col"), col.vert = par("col"), 
             col.main = par("col.main"), col.inner = par("col.main"), 
             col.sub = par("col.sub"),  cex.points = 2.0, 
             pch.u = 21, pch.a = 16, mfColRow = TRUE,
             to.draw.arg = NULL, withSubst = TRUE){

      xc <- match.call(call = sys.call(sys.parent(1)))$x
      ### manipulating the ... - argument
      dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."

      to.draw <- 1:3
      names(to.draw) <- c("d","p","q")
      if(! is.null(to.draw.arg)){
         if(is.character(to.draw.arg)) 
            to.draw <- pmatch(to.draw.arg, names(to.draw))
         else if(is.numeric(to.draw.arg)) 
              to.draw <- to.draw.arg
      }
      l.draw <- length(to.draw)

      pF <- expression({})
      if(!is.null(dots[["panel.first"]])){
         pF <- .panel.mingle(dots,"panel.first")
      }
      pF <- .fillList(pF, l.draw)
      pL <- expression({})
      if(!is.null(dots[["panel.last"]])){
          pL <- .panel.mingle(dots,"panel.last")
      }
      pL <- .fillList(pL, l.draw)
      dots$panel.first <- dots$panel.last <- NULL

      dots$ngrid <- NULL

      dots.for.points <- .makedotsPt(dots)
      dots.lowlevel <- .makedotsLowLevel(dots)
      dots.without.pch <- dots.lowlevel[! (names(dots.lowlevel) %in% c("col", "pch"))]
      ###
     if(!is(x,"DiscreteDistribution"))
         x <- .ULC.cast(x)     
      
     if(!is.logical(inner))
         {if(!is.list(inner))
              inner <- as.list(inner)
            #stop("Argument 'inner' must either be 'logical' or a 'list'")
          inner <- .fillList(inner,l.draw)          
         }

     cex <- if (hasArg("cex")) dots$cex else 1

     if (hasArg("cex") && missing(cex.points)) 
         cex.points <- 2.0 * cex

     if (hasArg("pch") && missing(pch.u))
          pch.u <- dots$pch
     if (hasArg("pch") && missing(pch.a))
          pch.a <- dots$pch
     
     if (hasArg("col") && missing(col.points))
         col.points <- dots$col
     if (hasArg("col") && missing(col.vert))
         col.vert <- dots$col
     if (hasArg("col") && missing(col.hor))
         col.hor <- dots$col
     if (hasArg("col") && missing(col.main))
        col.main <- dots$col
     if (hasArg("col") && missing(col.inner))
        col.inner <- dots$col
     if (hasArg("col") && missing(col.sub))
        col.sub <- dots$col

     if (!withSweave){
           devNew(width = width, height = height)
           }
     omar <- par("mar", no.readonly = TRUE)
 #    omar$cin <- omar$cra <- omar$csi <- omar$cxy <-  omar$din <- NULL
     if(mfColRow) on.exit(par(omar, no.readonly = TRUE))
     
     mainL <- FALSE
     subL <- FALSE
     lineT <- NA
     logpd <- logq <- ""

      ## getting the parameter
     slots <-  slotNames(param(x))
     slots <-  slots[slots != "name"]
     nrvalues <-  length(slots)
     if(nrvalues > 0){
           values <-  numeric(nrvalues)
       for(i in 1:nrvalues)
         values[i] <-  attributes(attributes(x)$param)[[slots[i]]]
       paramstring <-  paste(values, collapse = ", ")
       nparamstring <-  paste(slots, "=", values, collapse = ", ")
       qparamstring <- paste("(",paramstring,")",sep="")
     }
     else paramstring <- qparamstring <- nparamstring <- ""


     .mpresubs <- if(withSubst){
             function(inx) 
                    .presubs(inx, c("%C", "%D", "%N", "%P", "%Q", "%A"),
                          c(as.character(class(x)[1]), 
                            as.character(date()), 
                            nparamstring, 
                            paramstring, 
                            qparamstring,
                            as.character(deparse(xc))))
            }else function(inx) inx

     xlab0 <- list("d"="x", "p"="q", "q"="p")
     iL <- 1:length(to.draw)
     .mp2 <- function(dlb = dots$xlab, lb0 = list("d"="x", "p"="q", "q"="p")){
              dlb0 <- eval(dlb)
              if (!is.null(dlb)){
              .mp <- if(is.list(dlb0)) function(x,i){
                                if(is.call(x)) x <- eval(x)
                                if(length(i)==0) return(NULL)
                                i <- min(i)
                                if(is.character(x[[i]])){
                                   return(as.character(eval(.mpresubs(x[[i]]))))
                                }else{
                                res <- .mpresubs(x[[i]])
                                if(length(res)==0) return(NULL)
                                if(is.call(res)) res <- res[-1]
                                return(res)}
                                }else function(x,i){
                                  if(length(x)==1) return(x[1])
                                  res <- x[i]
                                  if(length(res)==0) return(NULL)
                                  if(is.na(res)) return(NULL)
                                  return(res)}
              force(lb0)
              .mp3 <- .mp(dlb,iL[to.draw==1])
              if(1%in%to.draw & !is.null(.mp3)) lb0[["d"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==2])
              if(2%in%to.draw & !is.null(.mp3)) lb0[["p"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==3])
              if(3%in%to.draw & !is.null(.mp3)) lb0[["q"]] <- .mp3

             }
             return(lb0)}
     xlab0 <- .mp2()
     dots$xlab <- NULL
     ylab0 <- .mp2(dlb=dots$ylab, lb0=list("d"="d(x)", "p"="p(q)", "q"="q(p)"))
     dots$ylab <- NULL

     if (hasArg("main")){
         mainL <- TRUE
         if (is.logical(main)){
             if (!main) mainL <-  FALSE 
             else 
                  main <- gettextf("Distribution Plot for %%A") ### 
                          ### double  %% as % is special for gettextf
             }
         main <- .mpresubs(main)
         if (mainL) {
             if(missing(tmar))
                tmar <- 5
             if(missing(cex.inner)) 
                cex.inner <- .9 
             lineT <- 0.6
             }
     }
     if (hasArg("sub")){ 
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
     
     if(mfColRow)
        opar <- par(mfrow = c(1,l.draw), mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE)
     else 
        opar <- par(mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE)
     
     if(is.logical(inner)){     
        inner.d <- if (inner) 
                   .mpresubs(gettextf("Probability function of %%C%%Q")) else 
                   ""
        inner.p <- if (inner) 
                   .mpresubs(gettextf("CDF of %%C%%Q")) else ""
        inner.q <- if (inner) 
                   .mpresubs(gettextf("Quantile function of %%C%%Q")) else ""
                          ### double  %% as % is special for gettextf
     }else{
        iL <- 1:length(to.draw)

        inner.d <- if(1%in%to.draw) .mpresubs(inner[[min(iL[to.draw==1])]]) else NULL
        inner.p <- if(2%in%to.draw) .mpresubs(inner[[min(iL[to.draw==2])]]) else NULL
        inner.q <- if(3%in%to.draw) .mpresubs(inner[[min(iL[to.draw==3])]]) else NULL
     }
                              
      lower <- min(support(x))
      upper <- max(support(x))
      dist <-  upper - lower

      supp <- support(x); 
      
      if(hasArg("xlim")) 
            {if(length(xlim) != 2) 
                stop("Wrong length of Argument xlim")
             supp <- supp[(supp >= xlim[1]) & (supp <= xlim[2])]         

            }else{
            }

       dx <- d(x)(supp)

     if(hasArg("ylim"))
         {if (3 %in% to.draw && any( c(1,2) %in% to.draw)){
                 if(! length(ylim) %in% c(2,4)) 
                     stop("Wrong length of Argument ylim")
           }else{
                 if(! length(ylim) == 2) 
                     stop("Wrong length of Argument ylim")
           }                  
           ylim <- matrix(ylim, 2,2)
           ylim1 <- ylim[,1]; ylim2 <- ylim[,2];
           }else{
              ylim1 <- c(0, max(dx)) 
              ylim2 <- c(-0.05,1.05)
              }

      if(hasArg("log"))
          {logpd <- dots$log
           logq <- gsub("u","y",gsub("y","x",gsub("x", "u", logpd)))
           if(length(grep("y",logpd))){ 
              ylim1 <- c(max(min(dx[dx>0]), ylim1[1]), 
                              ylim1[2])
              ylim2 <- c(max(min(dx[dx>0]), ylim2[1]), 
                              ylim2[2])
              }
           }

       o.warn <- getOption("warn")
       options(warn = -1)
       on.exit(options(warn=o.warn))

     plotCount <- 1

     if(1%in%to.draw){
       dots.without.pch$panel.first <- pF[[plotCount]]
       dots.without.pch$panel.last  <- pL[[plotCount]]
       do.call(plot, c(list(x = supp, dx, type = "h", pch = pch.a,
            ylim = ylim1, xlim=xlim, ylab = ylab0[["d"]], xlab = xlab0[["d"]],
            log = logpd), dots.without.pch))
       dots.without.pch$panel.first <- dots.without.pch$panel.last <- NULL
       plotCount <- plotCount + 1
       options(warn = o.warn)


       title(main = inner.d, line = lineT, cex.main = cex.inner,
             col.main = col.inner)

       if(do.points)
          do.call(points, c(list(x = supp, y = dx, pch = pch.a, 
                  cex = cex.points, col = col.points), dots.for.points))
       
       options(warn = -1)
       }
     ngrid <- length(supp)

     supp1 <- if(ngrid>1) supp else c(-max(1,abs(supp))*.08,0)+supp
     psupp1 <- c(0,p(x)(supp1))

     if(2%in%to.draw){
       dots.without.pch$panel.first <- pF[[plotCount]]
       dots.without.pch$panel.last  <- pL[[plotCount]]
       do.call(plot, c(list(x = stepfun(x = supp1, y = psupp1),
                     main = "", verticals = verticals, 
                     do.points = FALSE, 
                     ylim = ylim2, ylab = ylab0[["p"]], xlab = xlab0[["p"]],
                     col.hor = col.hor, col.vert = col.vert, 
                     log = logpd), dots.without.pch))
       dots.without.pch$panel.first <- dots.without.pch$panel.last <- NULL
       plotCount <- plotCount + 1
       if(do.points)
          {if(ngrid>1){
              do.call(points, c(list(x = supp, y = psupp1[1:ngrid], pch = pch.u, 
                  cex = cex.points, col = col.points), dots.for.points))
              do.call(points, c(list(x = supp, y = psupp1[2:(ngrid+1)], pch = pch.a, 
                  cex = cex.points, col = col.points), dots.for.points))
              }else{
              do.call(points, c(list(x = supp, y = 0, pch = pch.u, 
                  cex = cex.points, col = col.points), dots.for.points))           
              do.call(points, c(list(x = supp, y = 1, pch = pch.a, 
                  cex = cex.points, col = col.points), dots.for.points))           
              }
           }       
       options(warn = o.warn)

       
       title(main = inner.p, line = lineT, cex.main = cex.inner, 
             col.main = col.inner)

       if(do.points)
          do.call(points, c(list(x = supp, 
                  y = c(0,p(x)(supp[-length(supp)])), pch = pch.u, 
                  cex = cex.points, col = col.points), dots.for.points))
     }  

     if(3%in%to.draw){
       options(warn = -1)
       dots.without.pch$panel.first <- pF[[plotCount]]
       dots.without.pch$panel.last  <- pL[[plotCount]]
       do.call(plot, c(list(x = stepfun(c(0,p(x)(supp)),
                            c(NA,supp,NA), right = TRUE), 
            main = "", xlim = ylim2, ylim = c(min(supp),max(supp)),
            ylab = ylab0[["q"]], xlab = xlab0[["q"]],
            verticals = verticals, do.points = do.points, 
            cex.points = cex.points, pch = pch.a, 
            col.points = col.points,
            col.hor = col.hor, col.vert = col.vert, 
            log = logq), dots.without.pch))
       dots.without.pch$panel.first <- dots.without.pch$panel.last <- NULL
       plotCount <- plotCount + 1
       options(warn = o.warn)

      
       title(main = inner.q, line = lineT, cex.main = cex.inner,
             col.main = col.inner)

       dots.without.pch0 <- dots.without.pch
       dots.without.pch0$col <- NULL

       do.call(lines, c(list(x = c(0,p(x)(supp[1])), y = rep(supp[1],2),  
                  col = col.vert), dots.without.pch0))           

       if(do.points)
          {do.call(points, c(list(x = p(x)(supp[-length(supp)]),
                  y = supp[-1], pch = pch.u, cex = cex.points, 
                  col = col.points), dots.for.points))
           do.call(points, c(list(x = 0, y = supp[1], pch = pch.u, 
                  cex = cex.points, col = col.points), dots.for.points))}           
        
       if(verticals && ngrid>1)
          {dots.without.pch0 <- dots.without.pch
           dots.without.pch0$col <- NULL

           do.call(lines, c(list(x = rep(p(x)(supp[1]),2), y = c(supp[1],supp[2]),  
                  col = col.vert), dots.without.pch0))
          }
       }                      
       
     if (mainL)
           mtext(text = main, side = 3, cex = cex.main, adj = .5, 
                 outer = TRUE, padj = 1.4, col = col.main)                            
       
     if (subL)
           mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                 outer = TRUE, line = -1.6, col = col.sub)                            
   return(invisible())
   }
)

# -------- DistributionList   ---------- #

setMethod("plot", signature(x =  "DistrList", y = "missing"),
    function(x, ...){ 
        for(i in 1:length(x)){
            devNew()
            plot(x[[i]],...)
        }
        return(invisible())
    })
