`imgplot` <-
function (z=coralRG$R[,1], layout=coralRG$printer, crit1 = 0.05,
            crit2 = crit1, key.side=2,
            lohi.colors=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
              "#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),
            nacolor="#FFFF00", boxplot.side=1, split="quantiles")
{
  "block2matrix" <-
    function(z, sr=3, sc=2, gr=2, gc=2){
      ## Assumes that values in the vector z are  in row major
      ## order within blocks of dimension sr x sc, with blocks
      ## in row major order within a gr x gc array of grids.
      ## Elements in the vector that is returned are in row
      ## major order wrt the sr*gr x sc*gc matrix of values on
      ## the slide. (It is given the dimensions of a matrix.)
      xy <- array(z, dim=c(sc, sr, gc, gr))
      xy <- aperm(xy, c(1,3,2,4))
      dim(xy) <- c(sc*gc, gr*sr)
      xy}
  quantile.na <- function (z, ...)
    {
      tmp <- !(is.na(z) | is.infinite(z))
      quantile(z[tmp], ...)
    }
  length.na <- function (z, ...)
    {
      tmp <- !(is.na(z) | is.infinite(z))
      length(z[tmp], ...)
    }
  if(is.matrix(z))warning("z is a matrix, You probably want a column vector")
  bplot <- function(z, boxplot.side=1){
    xrange <- range(z,na.rm=TRUE)
    iqr <- diff(quantile(xrange, c(.25,.75)))
    bwex <- diff(xrange)/(3*iqr)
    xhi <- max(z,na.rm=TRUE)
    xusr <- par()$usr[c(1:2)]
    xpos=pretty(z[!is.na(z)], n=5)
    z <- xusr[1]+(z-xrange[1])*diff(xusr)/diff(xrange)
    newpos <- xusr[1]+(xpos-xrange[1])*diff(xusr)/diff(xrange)
    par(xpd=TRUE)
    atvert <- switch(boxplot.side, par()$usr[3]-par()$cxy[2]*0.8,
                     "", par()$usr[4]+par()$cxy[2]*0.8, "")
    if(atvert!=""){
      boxplot(z, at=atvert, boxwex=bwex, add=TRUE, horizontal=TRUE, xaxt="n")
      axis(side=boxplot.side, line=1.5,
           at=newpos, labels=xpos, cex.axis=0.75, mgp=c(2, 0.5, 0))
    }
    par(xpd=FALSE)
  }
  if (crit1 >= 1)
    crit1 <- crit1/(length.na(z))
  if (crit2 >= 1)
    crit2 <- crit2/(length.na(z))
  tmpind <- (z > quantile.na(z, probs = 1 - crit2)) | (z <
                                  quantile.na(z, probs = crit1))
  n <- prod(unlist(layout))
  n.all <- length(z)
  n.na <- sum(is.na(z))
  nhalf <- length(lohi.colors)%/%2
  n2 <- 2*nhalf
  n.one <- length(lohi.colors)
  plo <- crit1*(0:nhalf)/nhalf
  phi <- 1-crit2*(nhalf:0)/nhalf
  quiles1 <- quantile.na(z, plo)
  quiles2 <- quantile.na(z, phi)
  if(split=="intervals"){
    quiles1[2:nhalf] <- quiles1[1]+(quiles1[nhalf+1]-quiles1[1])*
      (1:(nhalf-1))/nhalf
    quiles2[2:nhalf] <- quiles2[1]-(quiles2[nhalf+1]-quiles2[1])*
      ((nhalf-1):1)/nhalf
    plo[-1] <- sapply(quiles1[-1],
                      function(x, z)sum(z<=x, na.rm=TRUE)/length.na(z), z=z)
    phi[-1] <- sapply(quiles2[-1],
                      function(x, z)sum(z<=x, na.rm=TRUE)/length.na(z), z=z)
  }

  if(crit1+crit2<1){
    quiles <- c(quiles1,quiles2)
    frac <- c(plo, phi)
    colpal <- c(lohi.colors[1:nhalf],"#FFFFFF",
                lohi.colors[(n.one-nhalf+1):(n.one)])
    midbreak <- TRUE
  }
  else {colpal <- lohi.colors
        midbreak <- FALSE
        quiles <- quantile.na(z, (0:n.one)/n.one)
        frac <- c(plo, phi[-1])
      }
  dups <- duplicated(quiles)
  if(any(dups)){
    cats <- seq(along=quiles[-1])
    filledcats <- cats[!dups]
    cutcats <- as.integer(cut(z, quiles[!dups], include.lowest=TRUE))
    fullm <- filledcats[cutcats]}
  else fullm <- as.integer(cut(z, quiles, include.lowest=TRUE))
  n.one <- length(colpal)
  nrects <- length(quiles)
  if(any(is.na(z))){
    nacat <- TRUE
    fullm[is.na(fullm)] <- max(unique(fullm[!is.na(fullm)]))+1
    colpal <- c(colpal, nacolor)
  }
  else nacat <- FALSE
  if ((length(as.vector(z)) != n) & (!is.null(names(z)))) {
    y <- fullm[tmpind]
    fullm <- rep(NA, n)
    fullm[as.integer(names(y))] <- y
  }
  else fullm[!tmpind] <- NA
  if ((length(as.vector(z)) != n) & (is.null(names(z)))) {
    stop(paste("Error: Length of vector is different from total number\n",
               "of spots and vector has no row.name.\n"))
  }
#################################################################
  gc <- layout$ngrid.c
  gr <- layout$ngrid.r
  sc <- layout$nspot.c
  sr <- layout$nspot.r
  full <- block2matrix(fullm, sr, sc, gr, gc)
  image(1:ncol(full), 1:nrow(full), t(full), axes = FALSE,
        xlab = "", ylab = "", col=colpal)
  box()
  abline(v = ((gr - 1):1) * (sr) + 0.5)
  abline(h = (1:(gc - 1)) * (sc) + 0.5)
#################################################################
  if(boxplot.side%in%c(1,3))bplot(z, boxplot.side=boxplot.side)
  if(key.side%in%c(2,4)){
    chw <- par()$cxy[1]
    barwid <- 0.75*chw
    if(key.side==2){
      x0 <- par()$usr[1]-chw-barwid
      xcutpos <- x0 - 0.4*chw
      xquilepos <- x0+barwid+0.55*chw
      srt <- 90
    }
    else {
      x0 <- par()$usr[2]+chw
      xcutpos <- x0 + barwid + 0.4*chw
      xquilepos <- x0-0.4*chw
      srt <- -90
    }
    yvals2 <- seq(from=par()$usr[3], to=par()$usr[4],
                  length=n2+midbreak+2*nacat+1)[-(n2+midbreak+2*nacat+1)]
    eps2 <- diff(yvals2[1:2])

    if(nacat){
      nlast <- length(yvals2)
      nclast <- length(colpal)
      rect(x0, yvals2[nlast], x0+barwid, yvals2[nlast]+eps2,
           col=colpal[nclast], xpd=TRUE)
      text(x0+0.5*barwid, yvals2[nlast]+0.5*eps2, "NA",
           xpd=TRUE, srt=srt)
      yvals2 <- yvals2[-((nlast-1):nlast)]
      colpal <- colpal[-nclast]
    }
    if(!midbreak){
      rect(x0, yvals2, x0+barwid, yvals2+eps2,
           col=colpal, xpd=TRUE)
      text(xcutpos, c(yvals2[1],yvals2+eps2),
           paste(signif(quiles,3)), srt=srt, xpd=TRUE, cex=0.8)
      text(xquilepos, yvals2[1], "(0%)", srt=srt, xpd=TRUE, cex=0.65)
      fracs <- frac[-c(1, length(frac))]
      text(xquilepos, yvals2[-1],
           paste("(",round(fracs*100,2),")",sep=""),
           srt=srt, xpd=TRUE, cex=0.65)
      text(xquilepos, yvals2[length(yvals2)]+eps2, "(100%)", srt=srt,
           xpd=TRUE, cex=0.65)

    }
    else {rect(x0, yvals2[1:nhalf], x0+barwid, yvals2[1:nhalf]+eps2,
               col=colpal[1:nhalf], xpd=TRUE)
          rect(x0, yvals2[(nhalf+2):(2*nhalf+1)], x0+barwid,
               yvals2[(nhalf+2):(2*nhalf+1)]+eps2,
               col=colpal[(nhalf+2):(2*nhalf+1)], xpd=TRUE)
          text(xcutpos, yvals2[1:(nhalf+1)],
               paste(signif(quiles1,3)), srt=srt, xpd=TRUE, cex=0.8)
          text(xquilepos, yvals2[2:(nhalf+1)],
               paste("(",round(plo[-1]*100,2),")",sep=""), srt=srt,
               xpd=TRUE, cex=0.65)
          text(xquilepos, yvals2[1], "(0%)", srt=srt, xpd=TRUE, cex=0.65)
          text(xcutpos,
               c(yvals2[(nhalf+2):(2*nhalf+1)], yvals2[2*nhalf+1]+eps2),
               paste(signif(quiles2,3)), srt=srt, xpd=TRUE, cex=0.8)
          text(xquilepos, yvals2[(nhalf+2):(2*nhalf+1)],
               paste("(",round(phi[-length(phi)]*100,2),")",sep=""),
               srt=srt, xpd=TRUE, cex=0.65)
          text(xquilepos, yvals2[2*nhalf+1]+eps2, "(100%)", srt=srt,
               xpd=TRUE, cex=0.65)
        }
  }
  invisible()
}

