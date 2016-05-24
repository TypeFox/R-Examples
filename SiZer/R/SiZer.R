SiZer <- function(x, y, h=NA, x.grid=NA, degree=NA, derv=1, grid.length=41){
  
  # calculate x.grid and h.grid from what was passed in.
  x.grid <- x.grid.create(x.grid, x, y, grid.length);
  h.grid <- h.grid.create(h, x, y, grid.length);
  
  # set up degree if it wasn't passed in
  if( is.na(degree) ){
  	 degree = derv+1;
  }
  
	row <- 1;
	slopes <- matrix(nrow=length(h.grid), ncol=length(x.grid));
	for( h in h.grid ){
		print(h)
		model <- locally.weighted.polynomial(x, y, h=h, x.grid=x.grid, degree=degree);     
  		intervals <- 
  	      calc.CI.LocallyWeightedPolynomial(model, derv=derv);
  		slopes[row,] <- find.states(intervals);         
		row <- row + 1;
	}

  out <- NULL;
  out$x.grid <- x.grid;
  out$h.grid <- h.grid;
  out$slopes <- slopes;
  class(out) <- 'SiZer';
  return(out);
}		



plot.SiZer <- function(x, ylab=expression(log[10](h)), 
			colorlist=c('red', 'purple', 'blue', 'grey'), ...){
	temp <- factor(x$slopes);
    final.colorlist <- NULL;
	if( is.element( '-1', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[1]);
	if( is.element( '0', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[2]);
	if( is.element( '1', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[3]);
	if( is.element( '2', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[4]);

	# Convert the slopes to a factor list to match up with final.colorlist
	temp <- matrix( as.integer(factor(x$slopes)), nrow=dim(x$slopes)[1] ) 

	image( x$x.grid, log(x$h.grid,10), t(temp), 
			col=final.colorlist, ylab=ylab, ...)
			
    # draw the bandwidth lines
    x.midpoint <- diff(range(x$x.grid))/2 + min(x$x.grid);
    lines( x.midpoint + x$h.grid, log(x$h.grid, 10), col='white' );
    lines( x.midpoint - x$h.grid, log(x$h.grid, 10), col='white' );                     
}	


h.grid.create <- function(h.grid, x, y, grid.length){
	foo <- grid.length;
	h.max <- diff( range(x) ) * 2;
	h.min <- max( diff(sort(x)) );
	if(all(is.na(h.grid))){
		# if we are passed an NA
		out <- 10^seq(log(h.min,10), log(h.max,10), length=grid.length);
	}else if( length(h.grid) == 1 ){
		# if we are passed a single integer
		out <- 10^seq(log(h.min,10), log(h.max,10), length=h.grid);
	}else if( length(h.grid) == 2 ){
		# passed two parameters: assume min and max for h.grid
		out <- 10 ^ seq(log(h.grid[1],10), 
                 		log(h.grid[2],10), length=grid.length);
	}else{
		# otherwise just return what we are sent
		out <- h.grid;
	} 
	return(out);
}


