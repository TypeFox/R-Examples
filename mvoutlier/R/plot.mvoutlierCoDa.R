plot.mvoutlierCoDa <-
function (x, ..., which=c("biplot","map","uni","parallel"),
          choice=1:2,coord=NULL,map=NULL,onlyout=TRUE,bw=FALSE,symb=TRUE,
          symbtxt=FALSE,col=NULL,pch=NULL,obj.cex=NULL,transp=1) 
{                                                          

# plot the object resulting from mvoutlier.CoDa
   if (onlyout){sel <- x$outliers} else {sel <- seq(1:length(x$outliers))}

    if (is.logical(symbtxt) & (!symbtxt)){
      typeask <- 1 # symbtxt is FALSE
      pchvec <- x$pchvec[sel]
    }
    else if (is.logical(symbtxt) & (symbtxt)){ # default text symbols
      typeask <- 2 # symbtxt is TRUE
      pchvec <- switch(onlyout+1,sel,c(1:sum(sel)))
    }
    else { # special text symbols provided in symbtxt
      typeask <- 2
      pchvec <- symbtxt[sel]
    }

# BIPLOT #############################
if (which[1]=="biplot"){

biplot.color <- function (x, y,mvobject=x,choicePC=choice,var.axes = TRUE, 
			  col, cex = rep(par("cex"), 2),
                          xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL,
                          ylim = NULL, arrow.len = 0.1, pch = NULL, main = NULL,
                          sub = NULL, axes = TRUE,
                          arrows.color = 1, box.color = 1,obj.cex=par("cex"), ...) {
  n <- nrow(x)
  p <- nrow(y)
  if (missing(xlabs)) {
    xlabs <- dimnames(x)[[1]]
    if (is.null(xlabs))
      xlabs <- 1:n
  }
  xlabs <- as.character(xlabs)
  dimnames(x) <- list(xlabs, dimnames(x)[[2]])
  if (missing(ylabs)) {
    ylabs <- dimnames(y)[[1]]
    if (is.null(ylabs))
      ylabs <- paste("Var", 1:p)
  }
  ylabs <- as.character(ylabs)
  dimnames(y) <- list(ylabs, dimnames(y)[[2]])
  if (length(cex) == 1)
    cex <- c(cex, cex)
  if (missing(col)) {
    col <- par("col")
    if (!is.numeric(col))
      col <- match(col, palette())
    box.color <- col
    arrows.color <- col + 1
  }
  else
    if (length(col) == 1) {
      box.color <- col
      arrows.color <- col
    }
  unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
  rangx1 <- unsigned.range(x[, 1])
  rangx2 <- unsigned.range(x[, 2])
  rangy1 <- unsigned.range(y[, 1])
  rangy2 <- unsigned.range(y[, 2])
  if (missing(xlim) && missing(ylim))
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
  else
    if (missing(xlim))
      xlim <- rangx1
  else ylim <- rangx2
  ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
  on.exit(par(oldpar))
  oldpar <- par(pty = "s")
  if (!is.null(main)) 
    oldpar <- c(oldpar, par(mar = par("mar") + c(0, 0, 1, 0)))
  if (missing(pch))
    temp.type <- "n"
  else
    temp.type <- "p"
# Plot-Symbol 
  eval <- mvobject$pcaobj$eigenvalues
  xperc <- eval[choicePC[1]]/sum(eval)*100
  yperc <- eval[choicePC[2]]/sum(eval)*100
  xaxlab <- paste("PC",choicePC[1]," (",round(xperc,1),"%)",sep="")
  yaxlab <- paste("PC",choicePC[2]," (",round(yperc,1),"%)",sep="")
  plot(x, type = temp.type, xlim = xlim, ylim = ylim,
       col = col, pch = pch,xlab=xaxlab,ylab=yaxlab,
       sub = sub, main = main, axes = axes,cex=obj.cex, ...)
# Text
  if (missing(pch))
    text(x, labels=xlabs, cex = cex[1], col = col, ...)
  par(new = TRUE)
  plot(y, axes = FALSE, type = "n", xlim = xlim * ratio,
       ylim = ylim * ratio, xlab = "", ylab = "", col = col,cex=obj.cex, ...)
  if (axes) {
    axis(3, col = arrows.color)
    axis(4, col = arrows.color)
  }
  box(col = box.color)
# Text  Variables
  text(y, labels = ylabs, cex = cex[2], col = arrows.color, ...)
  if (var.axes)
# Color of variable arrows
    arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8,
           col = arrows.color, length = arrow.len)
  invisible()
}

##############################

    pca <- x$pcaobj$princompOutputClr

    if (symb) {                                                
        if (bw){colv <- x$colbw} else {colv <- x$colcol}
          if (typeask==1){
            biplot.color(pca$scores[sel,choice],pca$loadings[,choice],
		mvobject=x,col=colv[sel], pch=x$pchvec[sel],obj.cex=x$cexvec[sel])
          }
          else{
            biplot.color(pca$scores[sel,choice],pca$loadings[,choice],
		mvobject=x,col=colv[sel],xlabs=pchvec,obj.cex=x$cexvec[sel])
          }
    }
    else {  # symb is FALSE 
        if (missing(pch)){pch <- (x$out+2)[sel]}
        if (missing(col)){
          if (bw){col <- ifelse(x$out[sel],gray(0.2),gray(0.7))} else {col <- (x$out+2)[sel]}
        }
        if (missing(obj.cex)){obj.cex <- 1}
          if (typeask==1){
            biplot.color(pca$scores[sel,choice],pca$loadings[,choice],
                mvobject=x,col=col, pch=pch,obj.cex=obj.cex)
          }
          else{
            biplot.color(pca$scores[sel,choice],pca$loadings[,choice],
                mvobject=x,col=col,xlabs=pchvec,obj.cex=obj.cex)
          }
    }
}

# PLOT MAP #############################
else if (which[1]=="map"){
    if (missing(coord)){stop("Please provide coordinates for the samples!")} 

    maprange <- apply(coord,2,range)
    if (symb) {
        if (bw){colv <- x$colbw} else {colv <- x$colcol}
          if (typeask==1){
            plot(coord[sel,],col=colv[sel], pch=x$pchvec[sel], cex=x$cexvec[sel],
                xlim=maprange[,1],ylim=maprange[,2],...) 
          }
          else{
            plot(coord[sel,],col=colv[sel], pch=x$pchvec[sel], cex=x$cexvec[sel],type="n",
                xlim=maprange[,1],ylim=maprange[,2],...) 
            text(coord[sel,],labels=pchvec,col=colv[sel],...) 
          }
    }
    else {  # symb is FALSE
        if (missing(pch)){pch <- (x$out+2)[sel]}
        if (missing(col)){col <- (x$out+2)[sel]}
        if (missing(obj.cex)){obj.cex <- 1}
          if (typeask==1){
            plot(coord[sel,],col=col, pch=pch, cex=obj.cex,
                xlim=maprange[,1],ylim=maprange[,2],...) 
          }
          else{
            plot(coord[sel,],col=col, pch=pch, cex=obj.cex,type="n",
                xlim=maprange[,1],ylim=maprange[,2],...) 
            text(coord[sel,],labels=pchvec,col=col,cex=obj.cex,...)
          }
    }
    if (!missing(map)){lines(map,col="gray")}
}
 
# Univariate Scatter Plot #############################
else if (which[1]=="uni"){
    Zj <- x$ilrvariables[sel,]
    ###old.par <- par(mfrow = c(1, ncol(Zj)), mai = c(0.6, 0, 0.6, 0), oma = c(0, 3, 0, 3))
    old.par <- par(mfrow = c(1, ncol(Zj)), mai = c(0.6, 0.15, 0.6, 0.15), oma = c(0, 3, 0, 3))
    if (symb) {
        if (bw){colv <- x$colbw} else {colv <- x$colcol}
          if (typeask==1){
            for (i in 1:ncol(Zj)) {                                          
              plot(runif(nrow(Zj), min = -1, max = 1), Zj[,i],
                main = paste("ilr(",dimnames(Zj)[[2]][i],")",sep=""), xlim = c(-1.5, 1.5),
		xlab = "",ylab="ilr-variables",xaxt = "n", col = colv[sel],
                pch=x$pchvec[sel], cex=x$cexvec[sel],...) 
              ###par(yaxt = "n")
            }
          } 
          else { # typeask is 2
            for (i in 1:ncol(Zj)) {
              plot(runif(nrow(Zj), min = -1, max = 1), Zj[,i],
                main = paste("ilr(",dimnames(Zj)[[2]][i],")",sep=""), xlim = c(-1.5, 1.5),
                xlab = "",ylab="ilr-variables",xaxt = "n", col = colv[sel],type="n",...)
              text(runif(nrow(Zj), min = -1, max = 1), Zj[,i], labels=pchvec,
                col = colv[sel],...)
              ###par(yaxt = "n")
            }
          }
    }
    else { # symb is FALSE
        if (missing(pch)){pch <- (x$out+2)[sel]}
        if (missing(col)){
          if (bw){col <- ifelse(x$out[sel],gray(0.2),gray(0.7))} else {col <- (x$out+2)[sel]}
        }
        if (missing(obj.cex)){obj.cex <- 1}
          if (typeask==1){
            for (i in 1:ncol(Zj)) {
              plot(runif(nrow(Zj), min = -1, max = 1), Zj[,i],
                main = paste("ilr(",dimnames(Zj)[[2]][i],")",sep=""), xlim = c(-1.5, 1.5),
                xlab = "",ylab="ilr-variables",xaxt = "n", col = col,
                pch=pch, cex=obj.cex,...)
              ###par(yaxt = "n")
            }
          }
          else { # typeask is 2
            for (i in 1:ncol(Zj)) {
              plot(runif(nrow(Zj), min = -1, max = 1), Zj[,i],
                main = paste("ilr(",dimnames(Zj)[[2]][i],")",sep=""), xlim = c(-1.5, 1.5),
                xlab = "",ylab="ilr-variables",xaxt = "n", col = col,type="n",...)
              text(runif(nrow(Zj), min = -1, max = 1), Zj[,i], labels=pchvec,
                col = col,cex=obj.cex,...)
              ###par(yaxt = "n")
            }
          }
    }
    par(old.par)
}


# Parallel coordinate plot #############################
else if (which[1]=="parallel"){
  # modified parcoord from MASS package:
  parcoordmy <- function (x, col = 1, lty = 1, var.label = FALSE, overplot = FALSE, subset=NULL, colsub=1, ltysub=1, ...) {
   # overplot=TRUE allows to again plot a subset of the data defined by the outliers
   rx <- apply(x, 2L, range, na.rm = TRUE)
   x <- apply(x, 2L, function(x) (x - min(x, na.rm = TRUE))/(max(x,
       na.rm = TRUE) - min(x, na.rm = TRUE)))
   matplot(1L:ncol(x), t(x), type = "l", col = col, lty = lty,
       xlab = "", axes = FALSE, ...)
   if (overplot){matlines(1L:ncol(x), t(x[subset,]),type = "l", col = colsub, lty = ltysub)}
   axis(1, at = 1L:ncol(x), labels = paste("ilr(",colnames(x),")",sep=""))
   for (i in 1L:ncol(x)) {
       lines(c(i, i), c(0, 1), col = "grey70")
       if (var.label)
           text(c(i, i), c(0, 1), labels = format(rx[, i], digits = 3),
               xpd = NA, offset = 0.3, pos = c(1, 3), cex = 0.7)
   }
   return(xnormed=x)
   }

    Zj <- x$ilrvariables[sel,]
    if (symb) {
      if (typeask==1){
        if (bw){colv <- x$colbw} 
        else { # adjust color with transparancy defined by transp
          colv <- x$colcol
          tr <- format(as.hexmode(round(transp*255)),width=2,upper.case=TRUE)
          colv[!is.na(colv)] <- paste(substr(colv[!is.na(colv)],1,7),tr,sep="")
        }
        ret <- parcoordmy(Zj,ylab="", col = colv[sel],...)
      }
      else { # typeask is 2
        if (bw){colv <- x$colbw; colvtr <- colv} 
        else { # adjust color with transparancy defined by transp
          colv <- x$colcol
          colvtr <- x$colcol
          tr <- format(as.hexmode(round(transp*255)),width=2,upper.case=TRUE)
          colvtr[!is.na(colv)] <- paste(substr(colv[!is.na(colv)],1,7),tr,sep="")
        }
        ret <- parcoordmy(Zj,ylab="", col = colvtr[sel],...)
        mtext(text=pchvec,side=2,line=runif(nrow(ret),-1.5,1.5),
              at=ret[,1],col = colv[sel],las=1,cex=0.7,...)
        mtext(text=pchvec,side=4,line=runif(nrow(ret),-1.5,1.5),
              at=ret[,ncol(ret)],col = colv[sel],las=1,cex=0.7,...)
      }
    }
    else { # symb is FALSE
      if (missing(pch)){pch <- (x$out+2)[sel]}
      if (missing(col)){
        if (bw){col <- ifelse(x$out[sel],gray(0.2),gray(0.7))} else {col <- (x$out+2)[sel]}
      }
      if (typeask==1){
        if (onlyout){
          ret <- parcoordmy(Zj,ylab="", col = col, ...)
        }
        else {
          # overplot outliers again:
          ret <- parcoordmy(Zj,ylab="", col = col, overplot=TRUE, subset=x$outliers, colsub=gray(0.2), ltysub=1, ...)
        }
      }
      else { # typeask is 2
        ret <- parcoordmy(Zj,ylab="", col = col,...)
        mtext(text=pchvec,side=2,line=runif(nrow(ret),-1.5,1.5),at=ret[,1],col = 1,las=1,cex=0.5,...)
        mtext(text=pchvec,side=4,line=runif(nrow(ret),-1.5,1.5),at=ret[,ncol(ret)],col = 1,las=1,cex=0.5,...)
      }
    }
}

else { stop("Plot not defined - select parameter WHICH accordingly!")}

invisible()

}

