get.uri.2d <-
function(x1, x2, tt, vv, spline.df=NULL){

  o <- order(x1, x2, decreasing=TRUE)
  
  # sort x2 by the order of x1
  x2.ordered <- x2[o]
  
  tv <- cbind(tt, vv)
  ntotal <- length(x1) # number of peaks    

  uri <- apply(tv, 1, comp.uri, x=x2.ordered)

  # compute the derivative of URI vs t using small bins
  uri.binned <- uri[seq(1, length(uri), by=4)]
  tt.binned <- tt[seq(1, length(uri), by=4)]
  uri.slope <- (uri.binned[2:(length(uri.binned))] - uri.binned[1:(length(uri.binned)-1)])/(tt.binned[2:(length(uri.binned))] - tt.binned[1:(length(tt.binned)-1)])

  # smooth uri using spline
  # first find where the jump is and don't fit the jump
  # this is the index on the left
  # jump.left.old  <- which.max(uri[-1]-uri[-length(uri)])
  short.list.length <- min(sum(x1>0)/length(x1), sum(x2>0)/length(x2))

  if(short.list.length < max(tt)){
    jump.left <- which(tt>short.list.length)[1]-1
  } else {
    jump.left <- which.max(tt)
  }

#  reversed.index <- seq(length(tt), 1, by=-1)
#  nequal <- sum(uri[reversed.index]== tt[reversed.index])
#  temp  <- which(uri[reversed.index]== tt[reversed.index])[nequal]
#  jump.left <- length(tt)-temp
 
  if(jump.left < 6){
   jump.left <- length(tt)
  }
    
 
  if(is.null(spline.df))
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=6.4)
  else{
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=spline.df)
  }
  # predict the first derivative
  uri.der <- predict(uri.spl, tt[1:jump.left], deriv=1)

  invisible(list(tv=tv, uri=uri, 
                 uri.slope=uri.slope, t.binned=tt.binned[2:length(uri.binned)], 
                 uri.spl=uri.spl, uri.der=uri.der, jump.left=jump.left,
                 ntotal=ntotal))
 }

