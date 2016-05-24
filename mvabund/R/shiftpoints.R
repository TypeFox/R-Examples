################################################################################
# function that returns a shift-vector or matrix with the same dimensions as x
# adding this shift to x or y makes it possible to see overlapping
# points (x_i,y_i) slightly shifted in a plot
################################################################################
shiftpoints <- function(x,y,sh=( max(x)-min(x))/100, centered=TRUE, method=1,
   reg=6, na.rm=TRUE)
{
    p<-NCOL(x)	   
    N<-NROW(x)
    
   if (p==1) {
      dims=N
      times=dims
   } else 
   {
      dims <- c(N,p)
      times=N*p
  }
  	
  x <- as.vector(as.matrix(x))
  y <- as.vector(as.matrix(y))
  if (length(y)!= length(x)) stop("x and y dimensions are different")
  z <- z.nona <- rep.int(0,times=times)
  nas <- NULL
  
	if(na.rm){
	x.nona<-attr(na.omit(x),"na.action")
	y.nona<-attr(na.omit(y),"na.action")
	nas <- unique(c(x.nona, y.nona))
	if(!is.null(nas)) {
		x <- x[-nas]
		y <- y[-nas]
		z.nona <- z[-nas]
		}
	} else if(any(is.na(x)) | any(is.na(y)) ) stop ("missing values")
  
  xuni<-unique(x)
  lxuni<-length(xuni)
  if (!lxuni==times)
  { 
   for (i in 1:lxuni)
    {
      xsame<-x==xuni[i]
      if (sum(xsame)>1)
      {
        ysa<-y[xsame]
        luniysa<-unique(ysa)
		luniysa<-luniysa
        for (j in 1: length(luniysa))
        {
          sameysa<- ysa==luniysa[j]
          if (sum(sameysa)>1)
           {
            # Get index of values in sameysa, shift each one different!
            xysame <- xsame
            xysame[xsame]<- sameysa
			
			if (method==1) {
        # method 1: if there are more than two:
        # shift is a fraction (relative to number of shifts) of sh
				div <- sum(sameysa)-1
		  }	else {
				# method 2: for up to reg, shift is sh/reg,
        # for more than reg: the same as method 1
				div <- max(reg-1,(sum(sameysa)-1))
      }
   
			z.nona[xysame]<-((1:sum(sameysa))-1)/div
			
			if (centered) z.nona[xysame] <- z.nona[xysame]-(sum(sameysa)-1)/(2*div)
			}
        }
      }
    }
  }
  if(!is.null(nas)) z[-nas] <- z.nona
   else z <- z.nona
  # Ensure that z has the same dim as x, y.
  z <- array(z*sh,dim=dims)
  
  return(z)
  
}

