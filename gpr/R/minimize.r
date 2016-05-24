minimize = function (X, f, .length, covfunc, x, y)
{
#Explanation from gpml Matlab package:
# SIG and RHO are the constants controlling the Wolfe-
# Powell conditions. SIG is the maximum allowed absolute ratio between
# previous and new slopes (derivatives in the search direction), thus setting
# SIG to low (positive) values forces higher precision in the line-searches.
# RHO is the minimum allowed fraction of the expected (from the slope at the
# initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
# Tuning of SIG (depending on the nature of the function to be optimized) may
# speed up the minimization; it is probably not worth playing much with RHO.

	toprint= FALSE
	INT = 0.1
	EXT = 3.0                
	MAX = 20                  
	RATIO = 10             
	SIG = 0.1
	RHO = SIG/2	 
	if(is.array( .length)){
		if (max(dim( .length)) == 2)
		{
			 red= .length[2] 
			 .length=  .length[1]
		}
	} else{ 
			red=1
	}
	if ( .length>0)
	{
		 S='Linesearch'
	} else{
		 S='Function evaluation'
	} 

	i.m= 0                                            # zero the run length counter
	ls_failed = 0                            # no previous line search has failed
	f_out =eval(call (f, X, covfunc , x, y, NULL , TRUE))
	f0= f_out[1][[1]]
	df0 = t(f_out[2][[1]]) #out put is a colum, R makes it a row when put it in list to output, I conver it back here
	fX = f0
	i.m = i.m + ( .length<0)                                            # count epochs?!
	s = -df0
	s=round(s*10000)/10000	
	d0 = -t(s)%*%s           # initial search direction (steepest) and slope
	x3 = red/(1-d0)                      # initial step is red/(|s|+1)
	mainloop = TRUE
	while( i.m < abs( .length) && mainloop)
	{                                      # while not finished
  		i.m = i.m + ( .length>0)                                     # count iterations?!
  		X0 = X
  		F0 = f0
  		dF0 = df0                  # make a copy of current values
  		if ( .length>0)
  		{
  			 	  M = MAX
  			} else{
  				  M = min(MAX, - .length- i.m)
  		}
  		whilerun = TRUE
 		while (whilerun ==TRUE)                             # keep extrapolating as long as necessary
    		{
    			x2 = 0
    			f2 = f0
    			d2 = d0
    			f3 = f0
    			df3 = df0
    			success = FALSE
    			while  (success == FALSE && M > 0)
    			{
    				M = M - 1
    				i.m = i.m + (.length<0)                         #count epochs?!
    				options(show.error.messages = FALSE)
    				 f_out2 =eval(call (f , X+ x3[1]*s, covfunc , x, y, NULL, TRUE))
        				f3=  f_out2[1][[ 1 ]][[ 1 ]]
        				df3 = t(f_out2[2][[1]])
        				f3=round(f3*10000)/10000
        				df3=round(df3*10000)/10000   #for sake of roanding to give answers near matlab
        				if (is.na(f3) || is.infinite(f3) || is.nan(f3)  || any(is.nan(df3) || is.na(df3) || is.infinite(df3))  )  
        				{
        					  cat( " ")
        					  x3 = (x2+x3)/2 
        				}else{
        					success = TRUE
        				  }
        			        	options(show.error.messages = TRUE)		  
        			}
        			if (f3 < F0){
        				X0 = X+x3[1]*s
        				F0 = f3
        				dF0 = df3
        			}         # keep best values
    			d3 = t(df3)%*%s                                       # new slope
    			if( d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0  )# are we done extrapolating?
     			{
     				 whilerun = FALSE
     				 break	 
    			}	
      			x1 = x2
      			f1 = f2
      			d1 = d2                       # move point 2 to point 1
    			x2 = x3 
    			f2 = f3
    			d2 = d3                        # move point 3 to point 2
    			A = 6*(f1-f2)+3*(d2+d1)*(x2-x1)                 # make cubic extrapolation
    			B = 3*(f2-f1)-(2*d1+d2)*(x2-x1)
    			x3 = x1-d1*(x2-x1)^2/(B+sqrt(abs(B*B-A*d1*(x2-x1)))) # num. error possible, ok!		
    			if ( (B*B-A*d1*(x2-x1)   < 0)[1]  || is.nan(x3) || is.infinite(x3) || x3 < 0) # num prob | wrong sign?
      			{	
      				x3 = x2*EXT                                 # extrapolate maximum amount
    			}else if( x3 > x2*EXT)                  # new point beyond extrapolation limit?
    			{
      				x3 = x2*EXT                                 # extrapolate maximum amount
      			}else if (x3 < x2+INT*(x2-x1) )         # new point too close to previous point?
      			{	
      				x3 = x2+INT*(x2-x1)
   			}
   			x3= round(x3*10000)/10000   #this is just to make same answers with matlab code
      		}
      		while ((abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0 ) # keep interpolating
      		{
    			if( d3 > 0 || f3 > f0+x3*RHO*d0)                         # choose subinterval
    			{
      				x4 = x3
      				f4 = f3
      				d4 = d3                      #move point 3 to point 4
      			}else{
      				x2 = x3
      				f2 = f3
      				d2 = d3                      # move point 3 to point 2
      			}
    			if (f4 > f0 ){          
      				x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2))  # quadratic interpolation
    			}else{
      				A = 6*(f2-f4)/(x4-x2)+3*(d4+d2)                    # cubic interpolation
      				B = 3*(f4-f2)-(2*d2+d4)*(x4-x2)
      				x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A        # num. error possible, ok!
    			}
      			
      			 if (is.nan(x3) || is.infinite(x3) )
      			 {
      				x3 = (x2+x4)/2               # if we had a numerical problem then bisect
    			 }
    			x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2))  # don't accept too close
    			f_out3 =eval(call (f , X+ x3*s, covfunc , x, y, NULL, TRUE))  
			f3 = f_out3[1] [[1]][[1]]
			df3 = t(f_out3[2][[1]])
    			if (f3 < F0)
    			{
    				x3=x3[[1]] # arrays of one in matlab converst to number, in R it is not
    				 X0 = X+x3*s
    				F0 = f3
    				dF0 = df3
    			}# keep best values
    			M = M - 1
    			i.m = i.m + (.length<0)                            # count epochs?!
    			d3 =  t(df3)%*%s                                                    # new slope
            }#end while                                        #end interpolation
 	   if (abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0 )        # if line search succeeded
 	   {
 	   	x3=x3[[1]]
    		X = X+x3*s
    		f0 = f3
    		fX=	t(cbind (t(fX), f0))
    		 s=(  (   (t(df3)%*%df3- t(df0)%*%(df3))[1]  )/ (  (t(df0)%*%(df0))[[1]]    ) *s  ) - df3
    		df0 = df3                                               # swap derivatives
    		d3 = d0
    		d0 = t(df0)%*%s
    		if (d0 > 0 ) 
    		{                                    # new slope must be negative
      			s = -df0
      			d0 = -t(s)%*%s                  # otherwise use steepest direction
    		}
    		x3 = x3 * min(RATIO, d3/(d0-  (  2^(-1022) )  ))          # slope ratio but max RATIO
   		 ls_failed = 0                           # this line search did not fail
 	  } else{
    		X = X0
    		f0 = F0
    		df0 = dF0                     # restore best point so far
   		 if (ls_failed || i.m > abs(.length)   )    # line search failed twice in a row
    		{
      			mainloop = 0 #break; 	                            #or we ran out of time, so we give up
      			break
    		}
    		s = -df0
     		d0 = -t(s)%*%s                                        #try steepest
    		x3 = 1/(1-d0)                    
    		ls_failed = 1                                    # this line search failed
  	  }#end if
	}#end first while
  return (list(X, fX, i.m))  #i.m is number of epochs
}