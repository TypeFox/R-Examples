trans.colors <-
function(col, alpha=0.5, names=NULL) {
  nc <- length(col)
  na <- length(alpha)
  # make lengths conform, filling out to the longest
  if (nc != na) {
  	col <- rep(col, length.out=max(nc,na))
  	alpha <- rep(alpha, length.out=max(nc,na))
  	}
  clr <-rbind(col2rgb(col)/255, alpha=alpha)
  col <- rgb(clr[1,], clr[2,], clr[3,], clr[4,], names=names)
  col
}

