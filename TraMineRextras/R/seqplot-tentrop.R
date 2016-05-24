## plot superposed transversal entropies
## author: Gilbert Ritschard

seqplot.tentrop <- function(seqdata, group,
     title=NULL, col=NULL, lty=NULL, lwd=3.5, ylim=NULL,
     xtlab=NULL, xtstep=NULL,
     withlegend=TRUE, glabels=NULL, legendpos="topright",
     horiz=FALSE, cex.legend=1, ...) {

  group <- factor(group)
  if(is.null(title)) {
     title <- "Tranversal Entropies"
  }

  entrop <- by(seqdata, group, seqstatd)

  k <- length(entrop)
  default.col <- brewer.pal(9,"Set1")
  ##default.col <- c("red","blue","black","magenta","green")
  if(is.null(col)) {
     #col <- colors.list[1:k]
     col <- default.col
  }
  kk <- ceiling(k/length(col))
  col <- rep(col,kk)
  col <- col[1:k]

  default.lty <- c("solid","dashed","dotted","solid","dashed")
  if(is.null(lty)) {
     lty <- default.lty
  }
  kk <- ceiling(k/length(lty))
  lty <- rep(lty,kk)
  lty <- lty[1:k]

  kk <- ceiling(k/length(lwd))
  lwd <- rep(lwd,kk)
  lwd <- lwd[1:k]

  npos <- ncol(seqdata)

  if(is.null(xtlab)){
      xtlab <- names(seqdata)
      ##xtlab <- paste("t",1:npos,sep="")
  }

  if(is.null(glabels)){
      glabels <- levels(group)
  }
  
  if(is.null(xtstep)){
      xtstep <- attr(seqdata,"xtstep")
  }
  
  
  if(is.null(ylim)){
        maxe <- max(entrop[[1]]$Entropy)
        mine <- min(entrop[[1]]$Entropy)
        for (i in 2:k){
            maxe <- max(maxe,entrop[[i]]$Entropy)
            mine <- min(mine,entrop[[i]]$Entropy)
        }
        ylim <- c(floor(10*mine),ceiling(10*maxe))/10
  }

  plot(0, type= "n", axes=FALSE, xlab="", ylab="Entropy", main=title, ylim=ylim, xlim=c(1,npos))
  for (i in 1:k) {
     lines(entrop[[i]]$Entropy, col=col[i],  type="l", lty=lty[i], lwd=lwd[i], ...)
  }
  axis(1,labels=xtlab[seq(from=1, to=npos, by=xtstep)],at=seq(from=1, to=npos, by=xtstep))
  axis(2)
  if(withlegend){
    legend(legendpos, legend=glabels,  lwd=lwd, lty=lty[1:k], col=col[1:k], horiz=horiz, cex=cex.legend)
  }

  return(k)
}


##k <- seqplot.tentrop(seqs.coh, group=seqs$highedu, xtlab=xtlab20)

##

seqplot.tentrop.m <- function(seqdata.list,
     title=NULL, col=NULL, lty=NULL, lwd=3.5, ylim=NULL,
     xtlab=NULL, xtstep=NULL,
     withlegend=TRUE, glabels=NULL, legendpos="topright",
     horiz=FALSE, cex.legend=1, ...) {

  ncurve <- length(seqdata.list)
  warn <- FALSE
  for (i in 2:ncurve){
    if(ncol(seqdata.list[[1]]) != ncol(seqdata.list[[i]])) {warn <- TRUE}
  }
  if (warn) {
    warning("sequence objects in seqdata.list are not all of same length")
    }
  if(is.null(title)) {
     title <- "Tranversal Entropies"
  }

  entrop <- lapply(seqdata.list, seqstatd)

  k <- length(entrop)
  default.col <- brewer.pal(9,"Set1")
  ##default.col <- c("red","blue","black","magenta","green")
  if(is.null(col)) {
     #col <- colors.list[1:k]
     col <- default.col
  }
  kk <- ceiling(k/length(col))
  col <- rep(col,kk)
  col <- col[1:k]

  default.lty <- c("solid","dashed","dotted","solid","dashed")
  if(is.null(lty)) {
     lty <- default.lty
  }
  kk <- ceiling(k/length(lty))
  lty <- rep(lty,kk)
  lty <- lty[1:k]

  kk <- ceiling(k/length(lwd))
  lwd <- rep(lwd,kk)
  lwd <- lwd[1:k]

  npos <- max(unlist(lapply(seqdata.list,ncol)))

  if(is.null(xtlab)){
      xtlab <- paste("t",1:npos,sep="")
  }

  if(is.null(glabels)){
      glabels <- paste("seq",1:length(seqdata.list),sep="")
  }
  
  if(is.null(xtstep)){
      xtstep <- attr(seqdata.list[[1]],"xtstep")
  }

  if(is.null(ylim)){
        maxe <- max(entrop[[1]]$Entropy)
        mine <- min(entrop[[1]]$Entropy)
        for (i in 2:k){
            maxe <- max(maxe,entrop[[i]]$Entropy)
            mine <- min(mine,entrop[[i]]$Entropy)
        }
        ylim <- c(floor(10*mine),ceiling(10*maxe))/10
  }

  plot(0, type= "n", axes=FALSE, xlab="", ylab="Entropy", main=title, ylim=ylim, xlim=c(1,npos))
  for (i in 1:k) {
     lines(entrop[[i]]$Entropy, col=col[i],  type="l", lty=lty[i], lwd=lwd[i], ...)
  }
  axis(1,labels=xtlab[seq(from=1, to=npos, by=xtstep)],at=seq(from=1, to=npos, by=xtstep))
  axis(2)
  if(withlegend){
    legend(legendpos, legend=glabels,  lwd=lwd, lty=lty[1:k], col=col[1:k], horiz=horiz, cex=cex.legend)
  }

  return(k)
}
