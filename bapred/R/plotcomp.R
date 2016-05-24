plotcomp <-
function(xp, groups, y, alpha, col=rainbow(length(groupsun)),
  xlim, ylim, ...) {

  if(missing(xlim)) xlim <- range(xp[,1])
  if(missing(ylim)) ylim <- range(xp[,2])

  groupsun <- levels(groups)
  if(missing(col)) {
    col <- rainbow(length(groupsun))
    colconv <- rainbow(length(groupsun), alpha=alpha)
  }
  else {
    colconv <- col
  }
  
  xpsub <- xp[groups==groupsun[1],]

  if(!missing(y)) {  
    yun <- levels(y)
    pchall <- as.character(as.numeric(as.numeric(y==yun[2])) + 1)
    plot(xpsub, xlim=xlim, ylim=ylim, col=col[1], cex=0.7, pch=pchall[groups==groupsun[1]], ...)
  }
  else
    plot(xpsub, xlim=xlim, ylim=ylim, col=col[1], pch=20, ...)


  dens2 <- MASS::kde2d(xpsub[,1], xpsub[,2], lims=c(min(xpsub[,1])-sd(xpsub[,1]), max(xpsub[,1])+sd(xpsub[,1]), min(xpsub[,2])-
    sd(xpsub[,2]), max(xpsub[,2])+sd(xpsub[,2])))
  contour(dens2, add=TRUE, col=colconv[1], nlevels=5)

  for(i in 2:length(groupsun)) {
    xpsub <- xp[groups==groupsun[i],]
	if(!missing(y))
      points(xpsub[,1], xpsub[,2], col=col[i], cex=0.7, pch=pchall[groups==groupsun[i]])
	else
      points(xpsub[,1], xpsub[,2], col=col[i], pch=20)
    dens2 <- MASS::kde2d(xpsub[,1], xpsub[,2], lims=c(min(xpsub[,1])-sd(xpsub[,1]), max(xpsub[,1])+sd(xpsub[,1]), min(xpsub[,2])-
      sd(xpsub[,2]), max(xpsub[,2])+sd(xpsub[,2])))
    contour(dens2, add=TRUE, col=colconv[i], nlevels=5)
  }

  return(list(col=col))

}
