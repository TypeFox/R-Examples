### ipairs.R

ipairs <- function(x,
                   pixs = 1,
                   zmax = NULL,
                   ztransf = function(x){x},
                   colramp = IDPcolorRamp,
                   cex = par("cex"),
                   lab.diag,
                   cex.diag = NULL,
                   main = NULL,
                   d.main = 1.5,
                   cex.main = 1.5*par("cex.main"),
                   legend = TRUE,
                   d.legend = 1.5,
                   cex.axis = 0.8*par("cex.axis"),
                   nlab.axis = 5,
                   minL.axis = 2,
                   las = 1,
                   border = FALSE,
                   mar = rep(0,4),
                   oma = c(3,3,1,0),
                   mgp = c(2,0.5,0)*cex.axis,
                   tcl = -0.3,
                   ...)
  ## Produces an image scatter plot matrix of largee datasets

  ## based on R function pairs V1.7
  ## Authors: Andreas Ruckstuhl, Rene Locher
  ## Version: 2009-04-08

{
  if (!(is.data.frame(x)||is.matrix(x)))
    stop("x must be a data.frame or a matrix\n")

  if(!all(sapply(x, function(x)
             any(is.element(is(x),c("numeric","factor","logical"))))))
    stop("All columns in data.frame must be either pure numerics, logicals or factors!\n")
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(oma=oma, mar=rep(0,4))

  nc <- ncol(x)
  if (missing(lab.diag)) {
    lab.diag <- colnames(x)
    if (is.null(lab.diag))
      lab.diag <- paste("var", 1:nc)
  }

  ## logicals are coerced into factors
  is.l <- sapply(x,is.logical)
  if (sum(is.l)>0) {
    for (i in which(is.l))
      x[,i] <- factor(x[,i],levels=c(FALSE,TRUE),labels=c("F","T"))
  }

  ## factors are sorted and coerced into integers
  is.f <- sapply(x,is.factor)
  nf <- sum(is.f)
  if(nf>0) {
    ##integer coercing of data.frames does not work!
    lev.list <- sapply(x[is.f],levels,simplify=FALSE)
    for (i in which(is.f)) x[,i] <- as.integer(x[,i])
    x <- cbind(x[,!is.f],x[,is.f])
    lab.diag <- c(lab.diag[!is.f],lab.diag[is.f])
    ## warning("factors have been put to the end and coerced into integers!\n")
  }

  is.f <- c(rep(FALSE,nc-nf),rep(TRUE,nf))
  w <- par("cin")[1] * 2.54
  h <- par("cin")[2] * 2.54

  ## 40% of space for color bar and 60% of it for axis labels
  w.legend <- lcm(7*cex.axis*cex*w)
  h.main <- lcm(cex.main*cex*h)
  d.main <- lcm(d.main*cex.main*cex*h)
  d.legend <- lcm(d.legend*cex.main*cex*h)

  if (!is.null(main) & legend) { ## plot title and legend
      lom <- matrix(1:nc^2, ncol=nc) + 2
      lom <- rbind(c(rep(1,nc),rep(0,2)), rep(0,nc+2),
                   cbind(lom,rep(0,nc),rep(2,nc)))
      lo <- layout(lom,
                   widths  = c(rep(1,nc),d.legend,w.legend),
                   heights = c(h.main,d.main,rep(1,nc)),
                   respect=TRUE)

      iplotMain(main,cex.main,cex=cex)
      iplotLegend(colramp=colramp,ncol=zmax,cex.axis=cex.axis,
                  border=border, mar=c(mar[1],0,mar[3],4)*cex.axis,
                  las=las, tcl=tcl, cex=cex)
  } ## plot title and legend

  if (is.null(main) & legend) { ## plot legend only
      lom <- matrix(1:nc^2, ncol=nc) + 1
      lom <- rbind(cbind(lom,rep(0,nc),rep(1,nc)))

      lo <- layout(lom,
                   widths  = c(rep(1,nc), d.legend, w.legend),
                   heights = rep(1,nc),
                   respect=TRUE)

       iplotLegend(colramp=colramp,ncol=zmax,cex.axis=cex.axis,
                  border=border, mar=c(mar[1],0,mar[3],4)*cex.axis,
                  las=las, tcl=tcl, cex=cex)
  }## plot legend only

  if (!is.null(main) & !legend) { ## plot title only
      lom <- matrix(1:nc^2, ncol=nc) + 1
      lom <- rbind(rep(1,nc), rep(0,nc),lom)

      lo <- layout(lom,
                   widths  = rep(1,nc),
                   heights = c(h.main, d.main, rep(1,nc)),
                   respect=TRUE)
      iplotMain(main,cex.main,cex=cex)
  } ## plot title only

  if (is.null(main) & !legend) { ## No title, no legend
      lom <- matrix(1:nc^2, ncol=nc)

      lo <- layout(lom,
                   widths  = rep(1,nc),
                   heights = rep(1,nc),
                   respect=TRUE)
  }  ## No title, no legend

##   layout.show(lo)
##   browser()

  cntsmax <- 0

  ## drawing scatterplots and axes
  par(mar=mar*cex.axis, las=las, cex=cex, cex.axis=cex.axis,
      mgp=mgp, pty="s", tcl=tcl, ...)
  for (i in 1:nc){
    for (j in 1:nc) {
      plot(if(is.f[i]) range(NaRV.omit(x[, i]))+0.5*c(-1,1) else
           range(NaRV.omit(x[, i])),
           if(is.f[j]) range(NaRV.omit(x[, j]))+0.5*c(-1,1) else
           range(NaRV.omit(x[, j])),
           xlab="", ylab="", axes=FALSE, type="n")
      box()

      if (i == 1 && (!(j%%2))) {## draw axes at top
        if(is.f[j]) {
          xmin <- min(x[, j],na.rm=TRUE)
          xmax <- max(x[, j],na.rm=TRUE)
          at <- seq(xmin, xmax,
                    by=max(floor((xmax-xmin)/(max(nlab.axis-1,1))),1))
          axis(2, at=at,
               labels=abbreviate(lev.list[[j-nc+nf]][at],
                                 minlength=minL.axis))
        } else
        axis(2, at=pretty(NaRV.omit(x[, j]),n=nlab.axis))
      }

      if (i == nc && (j%%2 )) {## draw axes at bottom
        if(is.f[j]) {
          xmin <- min(x[, j],na.rm=TRUE)
          xmax <- max(x[, j],na.rm=TRUE)
          at <- seq(xmin, xmax,
                    by=max(floor((xmax-xmin)/(max(nlab.axis-1,1))),1))
          axis(4, at=at,
               labels=abbreviate(lev.list[[j-nc+nf]][at],
                                 minlength=minL.axis))
        } else
        axis(4, at=pretty(NaRV.omit(x[, j]),n=nlab.axis), adj=1)
      }

      if (j == 1 && (!(i%%2))) {## draw axes at right side
        if(is.f[i]) {
          xmin <- min(x[, i],na.rm=TRUE)
          xmax <- max(x[, i],na.rm=TRUE)
          at <- seq(xmin, xmax,
                    by=max(floor((xmax-xmin)/(max(nlab.axis-1,1))),1))
          axis(3, at=at,
               labels=abbreviate(lev.list[[i-nc+nf]][at],
                                 minlength=minL.axis))
        } else
        axis(3, at=pretty(NaRV.omit(x[, i]),n=nlab.axis))
      }

      if (j == nc && (i%%2 ))  {## draw axes at left side
        if(is.f[i]) {
          xmin <- min(x[, i],na.rm=TRUE)
          xmax <- max(x[, i],na.rm=TRUE)
          at <- seq(xmin, xmax,
                    by=max(floor((xmax-xmin)/(max(nlab.axis-1,1))),1))
          axis(1, at=at,
               labels=abbreviate(lev.list[[i-nc+nf]][at],
                                 minlength=minL.axis))
        } else
        axis(1, at=pretty(NaRV.omit(x[, i]),n=nlab.axis))
      }

      if(i!=j){## do scatter plot
        cntsmax <-
          max(cntsmax,
              Image(x=x[, i], y=x[, j],
                    pixs=pixs, zmax=zmax, ztransf=ztransf,
                    colramp=colramp, factors=c(is.f[i],is.f[j])))
      }

      else{## fill text into diagonal
        par(usr = c(0, 1, 0, 1))
        if (is.null(cex.diag)) {
          cex.diag <- min(cex.main, 0.8/max(strwidth(lab.diag, "user")))
      }
        text(x=0.5, y=0.5, labels=lab.diag[i], cex = cex.diag, font = 1)
      }
    }
  }

  invisible(c(cntsmax=cntsmax))
} # ipairs

