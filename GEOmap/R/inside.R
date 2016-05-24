`inside` <-
function(A, POK )
  {
	FLTMAX = 1E+37
  #  use the point to test to build a line with one end outside of array */
	lt  = list(p1=A, p2=A)
	lt$p2$x = FLTMAX
	lp  = list(p1=list(x=POK$x, y=POK$y), p2=list(x=POK$x, y=POK$y))
	j = 1 
	count = 0
	for( i in 2:length(POK$x) )
 	{
	lp  = list(p1=list(x=POK$x[i], y=POK$y[i]), p2=list(x=POK$x[i], y=POK$y[i]))
	# print(paste(sep=" ", lp$p1$x, lp$p1$y))
   	 if(!Lintersect(lp, lt))
     	 {
     	 lp$p2  = list(x=POK$x[j], y=POK$y[j])
 # 	print(paste(sep=" ", i, j, lp$p1$x, lp$p1$y,lp$p2$x,lp$p2$y    ))
     	 j = i 
 #       /* check if last point not on line intersects */
     	 if(Lintersect(lp, lt)) count=count+1 
     	 }
	}
	K = count-2*floor(count/2)	
  return(K) 
}

