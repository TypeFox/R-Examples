plot.formLatticeOutput <-
function(x,...){
  #  
  #  This function plots the lattice.  If there are
  #  unwanted lattice links, you may edit the lattice using
  #  editFromLatticeOutput
  #
  require(spdep)
  if(class(x)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  nodes <- x$nodes
  poly <- x$poly
  latt <- x$latt
  xp.min <- min(poly[,1])
  xp.max <- max(poly[,1])
  yp.min <- min(poly[,2])
  yp.max <- max(poly[,2])
  plot(coords=nodes,latt,xlim=c(xp.min,xp.max), ylim= c(yp.min, yp.max),
    cex=0.2,...)
  lines(rbind(poly,poly[1,]),lwd=1,...)
  if(min(card(latt))==0){
      points(nodes[card(latt)==0,],pch=19,...)
  }
}

