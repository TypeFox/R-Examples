
plot.cem.plot <- function( ps.plot ){

  ps <- ps.plot$ps
  path <- ps.plot$path
  dev1 <- ps.plot$dev1
  dev2 <- ps.plot$dev2
  fy <- ps.plot$fy
  nsegs <- ps.plot$nsegs

  #compute path
  if(!is.null(path)){

    n <- nsegs * (nrow(path)-1) + 1
    px <- matrix(nrow=n, ncol=ncol(path))
    for(i in 2:nrow(path) ){
      s <- nsegs*(i-2)+1
      e <- nsegs*(i-1)+1
      for(j in 1:ncol(path)){
        px[s:e,j] <- seq(path[i-1,j], path[i,j], length.out=nsegs+1)
      }
    }


    py <- predict(ps, px) 
    pv <- predict(ps, px, type="variance")
    pdx <- predict(ps, px, type="densityX")
    pl <- c(0)
    el <- c(0)
    for(i in 2:nrow(py$y) ){
      pl[i] = pl[i-1] + sqrt( sum( (py$y[i, ] - py$y[i-1, ])^2 ) )
      el[i] <- sqrt( sum((py$y[1, ] - py$y[i, ])^2 ) )
    }
    xl <- seq(0, sqrt( sum( (px[1, ] - px[n, ])^2 ) ), length.out=n)

    #compute dists to path
    dp <- matrix(nrow=n, ncol=nrow(fy))
    for(i in 1:n){ 
      dp[i, ] = sqrt( colSums( apply(fy, 1, "-", px[i, ])^2 ) )
    }
    dpmin <- apply(dp, 2, min)
  }

  if( !is.null(path) ){
    dl <- dev.list()
    if( (length(dl) == 0 ) || is.null(dev2) ){
      dev.new()
      dev2 <- dev.cur()
    }
    else if(is.element(dev2, dev.list()) ){
      dev.set(dev2)
    }
    else{
      dev.new()
      dev2 <- dev.cur()
    }
    if(dev2 == dev1){
      dev1 = NULL
    }
    #plot layout
    layout(matrix(c(1,2,5,5,3,4,6,6), 4, 2))
    
    #Variance plot along path
    plot(xl, sqrt(pv), col = "gray", type="l", lwd=3,
         xlab=expression(x), ylab=expression(sigma[projection]), 
         cex.lab=1.5, cex.axis=1.5, bty="n")

    plot(xl, pdx, col = "gray", type="l", lwd=3,
               xlab=expression(x), ylab="density", ylim = c(0, max(pdx)), 
               cex.lab=1.5, cex.axis=1.5, bty="n")


    #Distances
    plot(xl, pl, col = "dodgerblue2", type="l", lwd=3,
         xlab=expression(x), ylab=expression(distances), 
         cex.lab=1.5, cex.axis=1.5, bty="n")
    lines(xl, el, col = "gold2", lwd=3)

    plot(el, pl, col = "dodgerblue2", type="l", lwd=3,
               xlab="euclidean", ylab="geodesic", 
               cex.lab=1.5, cex.axis=1.5, bty="n")
    lines(el, el, col = "darkgrey", lwd=2, lty=2)


    #Manifold path plot
    image(z=py$y, y=(0:ncol(py$y))+0.5, x=xl, col = heat.colors(n = 100),
	  xlab="expression", ylab="sample" ) 
  

    #Linear path plot
    li <- matrix(nrow=n, ncol=ncol(py$y))
    for(j in 2:nrow(path) ){
     for(i in 1:ncol(py$y) ){
       s <- (j-2)*nsegs + 1
       e <- (j-1)*nsegs + 1
       li[s:e,i] =  seq(py$y[s, i], py$y[e, i], length.out=nsegs+1)
     }
    }
    image(z = li, col = heat.colors(n = 100), x=xl, zlim = range(py$y)) 
 
  }



  if(is.null(fy)){
    fy <- predict(ps)
  }
  
  nc <- ncol(fy)

  dl <- dev.list()
  if( (length(dl) == 0 ) || is.null(dev1) ){
    dev.new();
    dev1 <- dev.cur()
  }
  else if(is.element(dev1, dev.list()) ){
    dev.set(dev1)
  }
  else{
    dev.new()
    dev1 <- dev.cur()
  }

  if(ncol(fy) == 1){
     one <- fy
     one[] <- 1
     plot(fy, one,  col = "#00000020",  
	     pch=19, asp=1, cex=1, bty="n", xaxt ="n", 
	     yaxt = "n", xlab="", ylab="")

  }
  else{

  close.screen(all.screens=T)
  split.screen(c(nc, nc))

  for(i in 1:nc){
    for(j in 1:nc){
      if(i == j){
        if(!is.null(path)){
	  screen((j-1)*nc + i)
          par(mar = c(1,1,1,1))
	  if(i == 1){
            plot(NULL, xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), xlab="", ylab="", yaxt="n", xaxt="n", bty="n")

	    a <- pl/el
	    draw.arc(x=0, y=0, radius=1, angle1=-a/2 - pi/2, angle2=a/2 - pi/2, col="gold2", lwd=3,n=1)
	    draw.arc(x=0, y=0, radius=1, angle1=-a/2 - pi/2, angle2=a/2 - pi/2, col="dodgerblue2", lwd=3)
	  }
	  else if(i==2){
            par( mar = c(2,2,2,2) )
	    pxm <- apply(dp, 1, min)
	    kr  <- loess(pxm ~ xl)
	    dpp <- predict(kr)
            plot(pl, dpp, col="dodgerblue2", type="l", lwd=2, xlab="xpath", ylab=expression(sigma[x]), bty="n")
	  }
	}
      }
      else{
        screen((j-1)*nc + i)
        par(mar = c(1,1,1,1))
        #par(bg = 'white')

        if(is.null(path) ){ 
          plot(fy[,i], fy[, j],  col = "#00000020",  
	     pch=19, asp=1, cex=1, bty="n", xaxt ="n", 
	     yaxt = "n", xlab="", ylab="")
	}
	else{
	  plot(fy[,i], fy[, j],  col = "#90909020",  
	     pch=19, asp=1, cex=1, bty="n", xaxt ="n", 
	     yaxt = "n", xlab="", ylab="")
 
          index3 <- which(dpmin < 3*ps.plot$ps$sigmaX)
          index1 <- which(dpmin < ps.plot$ps$sigmaX)
          index3 <- setdiff(index3, index1)

	  points(fy[index3,i], fy[index3, j], col="#FFa50020", pch=19)
	  points(fy[index1,i], fy[index1, j], col="#FF450030", pch=19)

	}
	if(i==1){
	  axis(2)
	}
	if(i==nc){
	  axis(4)
	}
	if(j==nc){
	  axis(1)
	}
	if(j==1){
	  axis(3)
	}
	      
        if(!is.null(path) ){
          lines(px[, i], px[, j], col="dodgerblue2", lwd=3)
          draw.circle(x=px[1,i], y=px[1,j], radius=ps$sigmaX, lwd=2, border="dodgerblue2")
    #points(px[seq(1, n, nsegs), ], pch=as.character(xl[seq(1, n, nsegs) ]) )
        
        }
      }
    }
  }

  }#end if(ncol(fy) != 1)

  ps.plot$fy <- fy
  ps.plot$dev1 <- dev1 
  ps.plot$dev2 <- dev2

  ps.plot
  
}

plot.cem <- function (ps, path=NULL, nsegs = 100, ...){
 
   ps.plot <- structure(list(ps=ps, path=path, nsegs=nsegs), class="cem.plot")
 
   ps.plot <- plot(ps.plot)

   ps.plot
}



cem.change.path <- function(ps.plot, i=1, j=2, x=j){
 dev.set(ps.plot$dev1)

 fy <- ps.plot$fy
 nc <- ncol(fy)
 path <- ps.plot$path

 screen((j-1)*nc + i, new=F)
 xy <- xy.coords(fy[,i], fy[, j])

 for(i in 1:nrow(path) ) {
   l <- locator(n=1)
   if(x == i){
     a <- l$x
   }
   else{
     a <- l$y
   }
   path[i,x] = a 
 }
 ps.plot$path = path
 ps.plot

}


cem.create.path <- function(ps.plot, n=2, i=1, j=2){

 dev.set(ps.plot$dev1)

 fy <- ps.plot$fy
 nc <- ncol(fy)
 
 if(nc > 1){
   screen((j-1)*nc + i, new=F)
   xy <- xy.coords(fy[,i], fy[, j])

   path <- matrix(nrow=n, ncol=nc)
   for(i in 1:n) {
    index <- identify(xy, n=1)
    path[i,] = fy[index, ]
   }
 }
 else{
   one <- fy
   one[] <- 1
   xy <- xy.coords(fy, one)

   path <- matrix(nrow=n, ncol=nc)
   for(i in 1:n) {
    index <- identify(xy, n=1)
    path[i] = fy[index]
   }
 }
 ps.plot$path = path
 ps.plot

}





