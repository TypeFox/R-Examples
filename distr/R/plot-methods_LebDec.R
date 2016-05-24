############################ plot #######################

setMethod("plot", signature(x = "AffLinUnivarLebDecDistribution", y = "missing"),
    function(x, width = 10, height = 5.5, withSweave = getdistrOption("withSweave"),
             xlim = NULL, ylim = NULL, ngrid = 1000,
             verticals = TRUE, do.points = TRUE,
             main = FALSE, inner = TRUE, sub = FALSE,
             bmar = par("mar")[1], tmar = par("mar")[3], ...,
             cex.main = par("cex.main"), cex.inner = 1.2,
             cex.sub = par("cex.sub"), col.points = par("col"),
             col.hor = par("col"), col.vert = par("col"),
             col.main = par("col.main"), col.inner = par("col.main"),
             col.sub = par("col.sub"),  cex.points = 2.0,
             pch.u = 21, pch.a = 16, mfColRow = TRUE, to.draw.arg = NULL,
             withSubst = TRUE){

      mc <- as.list(match.call(call = sys.call(sys.parent(1)), expand.dots = TRUE)[-1])
      do.call(getMethod("plot",
              signature(x="UnivarLebDecDistribution",y="missing")), args = mc)
      return(invisible())
})

setMethod("plot", signature(x = "UnivarLebDecDistribution", y = "missing"),
    function(x, width = 10, height = 14.5, withSweave = getdistrOption("withSweave"),
             xlim = NULL, ylim = NULL, ngrid = 1000,
             verticals = TRUE, do.points = TRUE,
             main = FALSE, inner = TRUE, sub = FALSE,
             bmar = par("mar")[1], tmar = par("mar")[3], ...,
             cex.main = par("cex.main"), cex.inner = 0.9,
             cex.sub = par("cex.sub"), col.points = par("col"),
             col.hor = par("col"), col.vert = par("col"),
             col.main = par("col.main"), col.inner = par("col.main"),
             col.sub = par("col.sub"),  cex.points = 2.0,
             pch.u = 21, pch.a = 16, mfColRow = TRUE, to.draw.arg = NULL,
             withSubst = TRUE){

      mc <- match.call(call = sys.call(sys.parent(1)), expand.dots = TRUE)[-1]
      xc <- mc$x

      ### manipulating the ... - argument
      dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."

      plotD <- getMethod("plot", signature(x = "DiscreteDistribution", 
                                           y = "missing"))
      plotC <- getMethod("plot", signature(x = "AbscontDistribution", 
                                           y = "missing"))

      to.draw <- 1:8
      names(to.draw) <- c("p","q","d.c","p.c","q.c","d.d","p.d","q.d")
      if(! is.null(to.draw.arg)){
         if(is.character(to.draw.arg)) 
            to.draw <- pmatch(to.draw.arg, names(to.draw))
         else if(is.numeric(to.draw.arg)) 
              to.draw <- to.draw.arg
      }
      l.draw <- length(to.draw)

      xlab0.d <- xlab0.c <- list("d"="x", "p"="q", "q"="p")
      ylab0.d <- ylab0.c <- list("d"="d(x)", "p"="p(q)", "q"="q(p)")

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

      plotCount <- 1
      if(!is(x, "UnivarLebDecDistribution"))
          x <- .ULC.cast(x)

      if(is(x,"DiscreteDistribution")){
         mcl <- as.list(mc)
         mcl$to.draw.arg <- (1:3)[( (6:8) %in%to.draw )] 
         mcl$ngrid <- NULL
         whichPFL <- mcl$to.draw.arg   
         mcl$panel.first <- pF[whichPFL]
         mcl$panel.last  <- pL[whichPFL]
         if(is.null(mcl$xlab)) mcl$xlab <- xlab0.d
         if(is.null(mcl$ylab)) mcl$ylab <- ylab0.d
         if(!is.logical(inner)){
                if(length(inner)!=3)
                   {inner <- .fillList(inner, 8)
                     mcl$inner <- inner[6:8]}
                }                          
         do.call(plotD, mcl)
         return(invisible())
      }
      
      if(is(x,"AbscontDistribution")){
         mcl <- as.list(mc)
         mcl$col.hor <- NULL
         if(is.null(mcl$xlab)) mcl$xlab <- xlab0.c
         if(is.null(mcl$ylab)) mcl$ylab <- ylab0.c
         mcl$to.draw.arg <- (1:3)[( (3:5) %in%to.draw )]
         whichPFL <- mcl$to.draw.arg   
         mcl$panel.first <- pF[whichPFL]
         mcl$panel.last  <- pL[whichPFL]
            if(!is.logical(inner)){
                if(length(inner)!=3)
                   {inner <- .fillList(inner, 8)
                     mcl$inner <- inner[6:8]}
                }                          
         do.call(plotC, as.list(mcl))
         return(invisible())
      }
      
      
      if(.isEqual(x@mixCoeff[1],0)){
         x <- x@mixDistr[[2]]
         mcl <- as.list(mc)
         if(is.null(mcl$xlab)) mcl$xlab <- xlab0.d
         if(is.null(mcl$ylab)) mcl$ylab <- ylab0.d
         mcl$x <- x
         mcl$to.draw.arg <- (1:3)[( (6:8) %in%to.draw )] 
         mcl$ngrid <- NULL
         whichPFL <- if(l.draw<8) mcl$to.draw.arg else 5+mcl$to.draw.arg  
         mcl$panel.first <- pF[whichPFL]
         mcl$panel.last  <- pL[whichPFL]
            if(!is.logical(inner)){
                if(length(inner)!=3)
                   {inner <- .fillList(inner, 8)
                     mcl$inner <- inner[6:8]}
                }                          
         do.call(plotD, as.list(mcl))
         return(invisible())
        }

      if(.isEqual(x@mixCoeff[1],1)){
         x <- x@mixDistr[[1]]
         mcl <- as.list(mc)
         if(is.null(mcl$xlab)) mcl$xlab <- xlab0.c
         if(is.null(mcl$ylab)) mcl$ylab <- ylab0.c
         mcl$x <- x
         mcl$to.draw.arg <- (1:3)[( (3:5) %in%to.draw )] 
         mcl$col.hor <- NULL
         whichPFL <- if(l.draw<8) mcl$to.draw.arg else 2+mcl$to.draw.arg  
         mcl$panel.first <- pF[whichPFL]
         mcl$panel.last  <- pL[whichPFL]
            if(!is.logical(inner)){
                if(length(inner)!=3)
                   {inner <- .fillList(inner, 8)
                     mcl$inner <- inner[6:8]}
                }                          
         do.call(plotC, as.list(mcl))
         return(invisible())
        }

      dots.for.points <- .makedotsPt(dots)

      dots.lowlevel <- .makedotsLowLevel(dots)
      dots.without.pch <- dots.lowlevel[! (names(dots.lowlevel) %in% c("pch", "col"))]
      dots.for.lines <- .makedotsL(dots)
      dots.v <- dots.for.lines
      dots.v$col <- NULL
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
 #    omar$cin <- omar$cra <- omar$csi <- omar$cxy <-  omar$din <- NULL
     if(mfColRow) (on.exit(par(omar, no.readonly = TRUE)))
     
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
                  }else function(inx)inx
 
     .mp2 <- function(dlb = dots$xlab, lb0 = list(list("p"="q", "q"="p"),
                          list("d"="x", "p"="q", "q"="p"),
                          list("d"="x", "p"="q", "q"="p"))){
              if (!is.null(dlb)){
              if(is.call(dlb)) dlb <- dlb[-1]
              .mp <- if(is.list(dlb)) function(x,i){
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
              iL <- 1:length(to.draw)
              force(lb0)
              .mp3 <- .mp(dlb,iL[to.draw==1])
              if(1%in%to.draw & !is.null(.mp3)) lb0[[1]][["p"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==2])
              if(2%in%to.draw & !is.null(.mp3)) lb0[[1]][["q"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==3])
              if(3%in%to.draw & !is.null(.mp3)) lb0[[2]][["d"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==4])
              if(4%in%to.draw & !is.null(.mp3)) lb0[[2]][["p"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==5])
              if(5%in%to.draw & !is.null(.mp3)) lb0[[2]][["q"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==6])
              if(6%in%to.draw & !is.null(.mp3)) lb0[[3]][["d"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==7])
              if(7%in%to.draw & !is.null(.mp3)) lb0[[3]][["p"]] <- .mp3
              .mp3 <- .mp(dlb,iL[to.draw==8])
              if(8%in%to.draw & !is.null(.mp3)) lb0[[3]][["q"]] <- .mp3
             }
             return(lb0)}

      xlab0 <- .mp2()
      xlab0.c <- xlab0[[2]]
      xlab0.d <- xlab0[[3]]
      dots$xlab <- NULL
      ylab0 <- .mp2(dlb = dots$ylab, lb0 = list(list("p"="p(q)", "q"="q(p)"),
                          list("d"="d(x)", "p"="p(q)", "q"="q(p)"),
                          list("d"="d(x)", "p"="p(q)", "q"="q(p)")))
      ylab0.c <- xlab0[[2]]
      ylab0.d <- ylab0[[3]]
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

     if(mfColRow){
        opar <- par("mfrow", mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE)
        ## common:
        com.drw <- (1:2)[( (1:2) %in%to.draw )]
        disc.drw <-(1:3)[( (3:5) %in%to.draw )]
        cont.drw <-(1:3)[( (6:8) %in%to.draw )]
        lcom.drw <- length(com.drw)
        ldisc.drw <- length(disc.drw)
        lcont.drw <- length(cont.drw) 
        
        nrw.drw <- (lcom.drw>0)+(ldisc.drw>0)+(lcont.drw>0)
        le.drw <- c(lcom.drw,ldisc.drw,lcont.drw)
        
        tw.drw <- any(le.drw==2); th.drw <- any(le.drw==3);
        cmm.drw <- (1+tw.drw)*(1+th.drw*2)
        ma.drw <- NULL; man.drw <- 0
        if(lcom.drw>0){
           ma.drw <- rep( 1:lcom.drw, 
                          each = cmm.drw/length(com.drw))
           man.drw <- length(com.drw)
           }
        if(ldisc.drw>0){
           ma.drw <- c(ma.drw, rep(1:ldisc.drw + man.drw, 
                       each = cmm.drw/length(disc.drw)))
           man.drw <- man.drw + length(disc.drw)
           }
        if(lcont.drw>0){
           ma.drw <- c(ma.drw, rep(1:lcont.drw + man.drw, 
                       each = cmm.drw/length(cont.drw)))
           }
        if(nrw.drw >0 ) layout(matrix(ma.drw, byrow=TRUE,  nrow=nrw.drw))
     }else 
        opar <- par(mar = c(bmar,omar[2],tmar,omar[4]), no.readonly = TRUE)

     if(is.logical(inner)){
        inner.p <- if (inner)
                   .mpresubs(gettextf("CDF of %%C%%Q")) else ""
        inner.q <- if (inner)
                   .mpresubs(gettextf("Quantile function of %%C%%Q")) else ""
                          ### double  %% as % is special for gettextf
     }else{
        iL <- 1:length(to.draw[to.draw<=2])

        inner.p <- if(1%in%to.draw) .mpresubs(inner[[min(iL[to.draw==1])]]) else NULL
        inner.q <- if(2%in%to.draw) .mpresubs(inner[[min(iL[to.draw==2])]]) else NULL
     }


     lower0 <- min(getLow(x@mixDistr[[1]],
                      eps = getdistrOption("TruncQuantile")*2),
                  getLow(x@mixDistr[[2]]))
     upper0 <- max(getUp(x@mixDistr[[1]],
                      eps = getdistrOption("TruncQuantile")*2),
                  getUp(x@mixDistr[[2]]))
     me <- q(x)(1/2); s <- q(x)(3/4)-q(x)(1/4)
     lower1 <- me - 6 * s
     upper1 <- me + 6 * s
     lower <- max(lower0, lower1)
     upper <- min(upper0, upper1)


     ## ngrid  nr of gridpoints
     ## exactq two p-values are considered equal if difference is
     ## is less than 10^-exactq in abs. value

     dist <- upper - lower
     del <- getdistrOption("DistrResolution")
     supp <- support(x)

     if(hasArg("xlim"))
     {  if(length(xlim)!=2) stop("Wrong length of Argument xlim");
           grid <- seq(xlim[1], xlim[2], length = ngrid)
           supp <- supp[(supp >= xlim[1]) & (supp <= xlim[2])]
     }else{grid <- seq(from = lower - 0.1 * dist, to = upper + 0.1 * dist,
                       length = ngrid)
           supp <- support(x)
     }

     grid <- unique(sort( c(supp, supp-del , grid )))
     pxg <- p(x)(grid)


     
     if(hasArg("ylim"))
         { if (any(c(2,5,8) %in% to.draw) && any( c(1,3,4,6,7) %in% to.draw)){
                 if(! length(ylim) %in% c(2,4)) 
                     stop("Wrong length of Argument ylim")
           }else{
                 if(! length(ylim) == 2) 
                     stop("Wrong length of Argument ylim")
           }                  
           ylim <- matrix(ylim, 2,2)
           ylim2 <- ylim[,2]
           }
     else ylim2 <- c(-0.05,1.05)

     if(hasArg("log"))
         {logpd <- dots$log
          logq <- gsub("u","y",gsub("y","x",gsub("x", "u", logpd)))
          if(length(grep("y",logpd))){
             ylim2 <- c(max(min(pxg[pxg>0]), ylim2[1]),
                               ylim2[2])
             }
          }

     if(!verticals){
         grid <- unique(sort( c(supp-del/2, grid )))
         grid[.isIn(grid,cbind(supp-del/2,supp-del/2))] <- NA
         pxg <- p(x)(grid)
     }else{
         xv <- as.vector(t(cbind(supp-del,supp,NA)))
         pxv <- p(x)(xv)
     }

     o.warn <- getOption("warn"); options(warn = -1)
     if(1 %in% to.draw){
        on.exit(options(warn=o.warn))
        dots.lowlevel$panel.first <- pF[[plotCount]]
        dots.lowlevel$panel.last  <- pL[[plotCount]]
        dots.lowlevel$xlim <- xlim
        do.call(plot, c(list(x = grid, pxg, type = "l",
             ylim = ylim2, ylab = ylab0[[1]][["p"]], xlab = xlab0[[1]][["p"]], log = logpd),
             dots.lowlevel))
        dots.lowlevel$panel.first <- dots.lowlevel$panel.last <- NULL
        dots.lowlevel$xlim <- NULL
        plotCount <- plotCount + 1
        options(warn = o.warn)
   
        pxg.d <- p(x)(supp)
        pxg.d0 <- p(x)(supp-del)
        if(do.points){
           do.call(points, c(list(x = supp, y = pxg.d, pch = pch.a,
                     cex = cex.points, col = col.points), dots.for.points))
           do.call(points, c(list(x = supp-del, y = pxg.d0, pch = pch.u,
                     cex = cex.points, col = col.points), dots.for.points))
        }
        if(verticals){
            do.call(lines, c(list(x = xv, y = pxv, col = col.vert),
                    dots.v))
        }
   
        title(main = inner.p, line = lineT, cex.main = cex.inner,
              col.main = col.inner)
     }
     ### quantiles

     ### fix finite support bounds
     ixg  <-  grid>=max(q(x)(0),lower) & grid <= min(q(x)(1),upper)
     pxg  <-   pxg[ixg]
     grid <-  grid[ixg]
     if(is.finite(q(x)(0))) {grid <- c(q(x)(0),grid); pxg <- c(0,pxg)}
     if(is.finite(q(x)(1))) {grid <- c(grid,q(x)(1)); pxg <- c(pxg,1)}

#     ### fix constancy regions of p(x)
#     if(isOldVersion(x)) x <- conv2NewVersion(x)
#
#     if(length(pxv))
#         do.call(lines, c(list(x = pxv, y = xv), dots.without.pch))
     
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

     if(2 %in% to.draw){
        options(warn = -1)
        dots.without.pch$panel.first <- pF[[plotCount]]
        dots.without.pch$panel.last  <- pL[[plotCount]]
        do.call(plot, c(list(x = po, xo, type = "n",
             xlim = ylim2, ylim = xlim, ylab = ylab0[[1]][["q"]], xlab = xlab0[[1]][["q"]],
             log = logq), dots.without.pch), envir = parent.frame(2))
        plotCount <- plotCount + 1
        dots.without.pch$panel.first <- dots.without.pch$panel.last <- NULL
        options(warn = o.warn)
   
   
        title(main = inner.q, line = lineT, cex.main = cex.inner,
              col.main = col.inner)
   
        options(warn = -1)
        do.call(lines, c(list(x=po, y=xo), dots.for.lines))
   #    if (verticals && !is.null(gaps(x))){
   #         do.call(lines, c(list(rep(pu1,2), c(gaps(x)[,1],gaps(x)[,2]),
   #                 col = col.vert), dots.without.pch))
   #     }
        options(warn = o.warn)
   
   
        if (verticals && !is.null(gaps(x))){
                pu <- rep(pu1,3)
                xu <- c(gaps(x)[,1],gaps(x)[,2],rep(NA,ndots))
                o <- order(pu)
                do.call(lines, c(list(pu[o], xu[o],
                        col = col.vert), dots.v))
         }
        if(!is.null(gaps(x)) && do.points){
            do.call(points, c(list(x = pu1, y = gaps(x)[,1], pch = pch.a,
                    cex = cex.points, col = col.points), dots.for.points) )
            do.call(points, c(list(x = pu1, y = gaps(x)[,2], pch = pch.u,
                    cex = cex.points, col = col.points), dots.for.points) )
       
        }

        if(do.points){
           if(is.finite(q(x)(0))) 
              do.call(points, c(list(x = 0, y = q(x)(0), pch = pch.u,
                   cex = cex.points, col = col.points), dots.for.points) )
           if(is.finite(q(x)(1))) 
              do.call(points, c(list(x = 1, y = q(x)(1), pch = pch.a,
                   cex = cex.points, col = col.points), dots.for.points) )
        }
        if (mainL)
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)
   
        if (subL)
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)
                  
        }
    mc.ac <- mc
    if(!is.logical(inner)) 
       mc.ac$inner <- lapply(inner[3:5], function(x) 
                             if(is.character(x))
                                as.character(eval(.mpresubs(x)))
                             else .mpresubs(x)) 
     mc.ac$xlab <- xlab0.c
     mc.ac$ylab <- ylab0.c
     mc.ac$mfColRow <- FALSE
     mc.ac$main <- FALSE
     mc.ac$sub <- FALSE
     mc.ac$x <- NULL 
     mc.ac$withSweave <- TRUE 
     mc.ac$to.draw.arg <- (1:3)[( (3:5) %in%to.draw )] 
     if(is.null(mc.ac$cex.inner))  mc.ac$cex.inner <- 0.9

     whichPFL <- plotCount-1+mc.ac$to.draw.arg  
     mc.ac$panel.first <- pF[whichPFL]
     mc.ac$panel.last  <- pL[whichPFL]

     do.call(plotC, c(list(acPart(x)),mc.ac), envir = parent.frame(2))
     plotCount <- plotCount + 3

     mc.di <- mc
     if(!is.logical(inner)) 
         mc.di$inner <- lapply(inner[6:8], function(x) 
                               if(is.character(x))
                                  as.character(eval(.mpresubs(x)))
                               else .mpresubs(x)) 
     mc.di$xlab <- xlab0.d
     mc.di$ylab <- ylab0.d
     mc.di$mfColRow <- FALSE
     mc.di$main <- FALSE
     mc.di$sub <- FALSE
     mc.di$x <- NULL
     mc.di$ngrid <- NULL
     mc.di$withSweave <- TRUE 
     mc.di$to.draw.arg <- (1:3)[( (6:8) %in%to.draw )]
     if(is.null(mc.di$cex.inner))  mc.di$cex.inner <- 0.9

     whichPFL <- plotCount-1+mc.di$to.draw.arg  
     mc.di$panel.first <- pF[whichPFL]
     mc.di$panel.last  <- pL[whichPFL]
     do.call(plotD, c(list(discretePart(x)),mc.di), envir = parent.frame(2))
     plotCount <- plotCount + 3
     return(invisible())
     
   }
   )

setMethod("plot", signature(x = "CompoundDistribution", y = "missing"),
           function(x, ...) {
           mc <- as.list(match.call(call = sys.call(sys.parent(1)), 
                            expand.dots = TRUE)[-1])
           do.call(getMethod("plot",signature(x = "UnivarLebDecDistribution", 
                                      y = "missing")),args=mc)
           return(invisible())
           })

