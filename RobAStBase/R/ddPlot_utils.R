.ddPlot.MatNtNtCoCo <- function(data, ...,
                                dist.x = NormType(), 
                                dist.y  = NormType(),
                                cutoff.x = cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x),
                                cutoff.y = cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y),
                                cutoff.quantile.x = 0.95,  
                                cutoff.quantile.y = cutoff.quantile.x,
                                transform.x,
                                transform.y = transform.x,
                                id.n,
                                lab.pts,
                                adj =0,
                                cex.idn = 1,
                                col.idn = par("col"),
                                lty.cutoff,
                                lwd.cutoff,
                                col.cutoff = "red",
                                text.abline = TRUE,
                                text.abline.x = NULL, text.abline.y = NULL,
                                cex.abline = par("cex"), col.abline = col.cutoff,
                                font.abline = par("font"), adj.abline = c(0,0),
                                text.abline.x.x = NULL, text.abline.x.y = NULL,
                                text.abline.y.x = NULL, text.abline.y.y = NULL,
                                text.abline.x.fmt.cx = "%7.2f",
                                text.abline.x.fmt.qx = "%4.2f%%",
                                text.abline.y.fmt.cy = "%7.2f",
                                text.abline.y.fmt.qy = "%4.2f%%",
                                jitt.fac = 10){

       dots <- match.call(expand.dots = FALSE)$"..."

       id.n1 <- 1:ncol(data)

       if(missing(id.n) || is.null(id.n))
          id.n <- id.n1


       if(missing(lab.pts)|| is.null(lab.pts)){
          lab.pts <-  if(!is.null(colnames(data))) colnames(data) else id.n1
       }

       data <- data[,id.n, drop = FALSE]
       lab.pts <- lab.pts[id.n]
       id.n1 <- id.n1[id.n]

       id.n0 <- 1:length(id.n)

       data.x <- if(missing(transform.x)||is.null(transform.x)) data else
                    transform.x(data)
                    #if(!is(IC.x,"IC")) stop("Argument 'IC.x' of 'ddPlot' must be an 'IC'")
                    #   apply(data,2,function(xx) evalIC(IC.x,xx))}
       data.y <- if(missing(transform.y)||is.null(transform.y)) data else
                    transform.y(data)
                    #if(!is(IC.y,"IC")) stop("Argument 'IC.y' of 'ddPlot' must be an 'IC'")
                    #   apply(data,2,function(xx) evalIC(IC.y,xx))}


      if(is.null(dist.x)) dist.x <- NormType()
      if(is.null(dist.y)) dist.y <- NormType()

      if(is.null(dots$xlab)) dots$xlab <- name(dist.x)
      if(is.null(dots$ylab)) dots$ylab <- name(dist.y)

      if(!is.null(dots$log)){
	          if(grepl("x",dots$log)) dots$xlab <- paste(dots$xlab, "(log-scale)",
	                                               sep="  ")
	          if(grepl("y",dots$log)) dots$ylab <- paste(dots$ylab, "(log-scale)",
	                                               sep="  ")
	       }

      if(is.null(cutoff.quantile.x))
         cutoff.quantile.x <- 0.95

      if(is.null(cutoff.quantile.y))
         cutoff.quantile.y <- cutoff.quantile.x

      if(is.null(cutoff.x))
         cutoff.x <- cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x)
      else {assign("norm", dist.x, environment(fct(cutoff.x)))
            assign("cutoff.quantile", cutoff.quantile.x, environment(fct(cutoff.x)))}

      if(is.null(cutoff.y))
         cutoff.y <- cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y)
      else {assign("norm", dist.y, environment(fct(cutoff.y)))
            assign("cutoff.quantile", cutoff.quantile.y, environment(fct(cutoff.y)))}

      if(!is(dist.x, "NormType")) stop("Argument 'dist.x' of 'ddPlot' must be of class 'NormType'")
      if(!is(dist.y, "NormType")) stop("Argument 'dist.y' of 'ddPlot' must be of class 'NormType'")
      if(!is(cutoff.x, "cutoff")) stop("Argument 'cutoff.x' of 'ddPlot' must be of class 'cutoff'")
      if(!is(cutoff.y, "cutoff")) stop("Argument 'cutoff.y' of 'ddPlot' must be of class 'cutoff'")

      ndata.x <- fct(dist.x)(data.x)
      ndata.y <- fct(dist.y)(data.y)
      
#      print(head(ndata.x))

      co.x <- fct(cutoff.x)(data.x)
      co.y <- fct(cutoff.y)(data.y)

      if(is.null(adj)) adj <- 0
      if(missing(cex.idn)||is.null(cex.idn))
         cex.idn <- if(is.null(dots$cex)) 1 else dots$cex

      if(missing(col.idn)||is.null(col.idn))
         col.idn <- if(is.null(dots$col)) par("col") else dots$col

      if(is.null(dots$lwd)) dots$lwd <- par("lwd")
      if(is.null(dots$lty)) dots$lty <- par("lty")

      col.cutoff <- rep(col.cutoff,length.out=2)
      if(missing(lty.cutoff) && !is.null(dots$lty)) lty.cutoff <- dots$lty
      if(missing(lwd.cutoff) && !is.null(dots$lwd)) lwd.cutoff <- dots$lwd
      if(missing(cex.abline) && !is.null(dots$cex)) cex.abline <- dots$cex
      if(missing(adj.abline) && !is.null(dots$adj)) lty.abline <- dots$adj
      if(missing(font.abline) && !is.null(dots$font)) font.abline <- dots$font

      pdots <- .makedotsLowLevel(dots)
      pdots$xlab <- dots$xlab
      pdots$ylab <- dots$ylab
      pdots$nsim <- NULL
      pdots$x <- NULL
      pdots$y <- NULL
      pdots$offset <- NULL
      pdots$pos <- NULL
      pdots$untf <- NULL

      abdots <- .makedotsAB(dots)
      if(!missing(lty.cutoff)) abdots$lty <- lty.cutoff[[1]]
      if(!missing(lwd.cutoff)) abdots$lwd <- lwd.cutoff[1]
      abdots$col <- col.cutoff[1]
      abdots$jitt.fac <- dots$jitt.fac

      abdots <- list(abdots,abdots)
      abdots$jitt.fac <- pdots$jitt.fac

      if(!is.null(abdots$lty))
	          if(is.list(lty.cutoff)) abdots[[2]]$lty <-  lty.cutoff[[2]]
      if(!is.null(abdots$lwd))
	         if(length(lwd.cutoff)>1) abdots[[2]]$lwd <-  lwd.cutoff[2]

      ab.textL <- rep(text.abline,length.out=2)
	    abtdots.x <- abtdots.y <- vector("list",0)
	    cex.abline <- rep(cex.abline, length.out = 2)
	    col.abline <- rep(if(!is.null(col.abline))
                          col.abline else "red", length.out = 2)
      font.abline <- rep(font.abline, length.out = 2)
      adj.abline <- matrix(rep(adj.abline,length.out=4),2,2)


	    .mpresubs <- function(inx)
                    .presubs(inx, c("%qx", "%qy", "%cx", "%cy"),
                        c(gettextf(text.abline.x.fmt.qx,
                             round(cutoff.quantile.x*100,1)),
                          gettextf(text.abline.y.fmt.qy,
                             round(cutoff.quantile.y*100,1)),
                          gettextf(text.abline.x.fmt.cx,
                             round(co.x,2)),
                          gettextf(text.abline.y.fmt.cy,
                          round(co.y,2))))
      
      if(!missing(lwd.cutoff)) abdots$lwd <- lwd.cutoff
      if(!missing(lty.cutoff)) abdots$lty <- lty.cutoff
      abdots$jitt.fac <- dots$jitt.fac

      abtdots.x$labels <- if(! is.null(text.abline.x))
                       .mpresubs(text.abline.x) else gettextf(
                              paste(text.abline.x.fmt.qx,"-cutoff = ",
	                                            text.abline.x.fmt.cx,sep=""),
                              cutoff.quantile.x*100,round(co.x,digits=2))
      abtdots.x$cex <- cex.abline[1]
	    abtdots.x$col <- col.abline[1]
	    abtdots.x$font <- font.abline[1]
	    abtdots.x$srt <- NULL
	    abtdots.x$adj <- adj.abline[,1]

      abtdots.y$labels <- if(! is.null(text.abline.y))
                       .mpresubs(text.abline.y) else gettextf(
                             paste(text.abline.y.fmt.qy,"-cutoff = ",
	                                            text.abline.y.fmt.cy,sep=""),
                             cutoff.quantile.y*100,round(co.y,digits=2))
	    abtdots.y$cex <- cex.abline[2]
	    abtdots.y$col <- col.abline[2]
	    abtdots.y$font <- font.abline[2]
	    abtdots.y$srt <- NULL
	    abtdots.y$adj <- adj.abline[,2]

      tdots <- .makedotsT(dots)
      tdots$cex <- cex.idn
      tdots$col <- col.idn
      tdots$offset <- dots$offset
      tdots$pos <- dots$pos
      tdots$adj <- adj

      pdots$log <- dots$log
      pdots$adj <- par("adj")

      adots <- pdots
      adots$col <- pdots$col.axis
      adots$lty <- pdots$lty.axis
      adots$adj <- par("adj")

      pdots$axes <- FALSE
      pdots$adj <- par("adj")
      ####

#      print(quantile(ndata.x))
#      print(co.x)
#      print(fct(cutoff.x))
#      print(dist.x)
#      print(get("norm", environment(fct(cutoff.x))))
#      print(quantile(ndata.y))
#      print(co.y)
#      print(fct(cutoff.y))
#      print(get("norm", environment(fct(cutoff.y))))
##

      if(!is.null(dots$xlim))
          if(is.logical(dots$xlim))
             if(dots$xlim) pdots$xlim <- c(min(ndata.x)/1.03,max(ndata.x,co.x)*1.03)
       if(!is.null(dots$ylim))
          if(is.logical(dots$ylim))
             if(dots$ylim) pdots$ylim <- c(min(ndata.y)/1.03,max(ndata.y,co.y)*1.03)


      id.x <- id.n0[ndata.x >= co.x*.999]
      id.y <- id.n0[ndata.y >= co.y*.999]
      id.xy <- intersect(id.x,id.y)


      id0.xy <- id.n1[id.xy]
      id0.x <- id.n1[id.x]
      id0.y <- id.n1[id.y]
      do.call(plot, args = c(list(x = ndata.x, y=ndata.y, type = "p"), pdots))
      do.call(box,args=c(adots))

      pusr <- par("usr")
      mid.x <- mean(pusr[c(1,2)])
      mid.y <- mean(pusr[c(3,4)])
      abtdots.y$x <- if(is.null(text.abline.y.x)) mid.x else text.abline.y.x
      abtdots.x$y <- if(is.null(text.abline.x.y)) mid.y else text.abline.x.y

      do.call(abline, args = c(list(v=co.x), abdots[[1]]))
	    do.call(abline, args = c(list(h=co.y), abdots[[2]]))

      if(ab.textL[1])
         do.call(text, args = c(list(y=co.y*1.03), abtdots.y))
#         do.call(text, args = c(list(co.x-5,mid.y,paste(cutoff.quantile.y*100,"%-cutoff = ",round(co.x,digits=2)),srt=90)))
      if(ab.textL[2])
         do.call(text, args = c(list(x=co.x*1.03), abtdots.x,srt=90))
#      do.call(text, args = c(list(mid.x,co.y+5,paste(cutoff.quantile.x*100," %-cutoff = ",round(co.y,digits=2)))))

      if(length(id.xy))
         do.call(text, args = c(list(jitter(ndata.x[id.xy],factor=jitt.fac),
                                     jitter(ndata.y[id.xy],factor=jitt.fac),
                                labels=lab.pts[id.xy]), tdots))
          #axis(side=4)
#      axis(side=1)

      return(list(id.x=id0.x, id.y= id0.y, id.xy = id0.xy,
             qtx = quantile(ndata.x), qty = quantile(ndata.y),
             cutoff.x.v = co.x, cutoff.y.v = co.y
             ))
}
