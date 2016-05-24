fourier <-
function( f, order=3, ... ) {
  # Fourier approximation is determined and graphed over the interval from 0 to 2*pi.
  # 'f':  The function to be approximated by Fourier analysis.
  # 'order' is the order of the Fourier transformation.
  # The numerical output consists of a0/2, a1, ..., an, b1, ..., b2.
  # Note: The equation is (constant) + a_1 cos x + ... + a_n cos(n x) + b_1 sin(x) + ... + b_n sin(n x)
  # $constant is constant term.
  # $cosine.coefficients are the cosine coefficients.
  # $sine.coefficients are the sine coefficients.
  # Examples: 
  #     >  fourier( function(x){ exp(-x)*(x-pi) }, 4 )
  #     >  fourier( function(x){ exp(-x) }, 7 )
  #     >  fourier( function(x){ (x-pi) }, 5 )
   if (!is.function(f))  stop("'f' must be a function.")
   if (!is.function(f))  stop("'f' must be a function.")
   if (!is.numeric(order))  stop("'order' must be an integer.")
   if (length(order)!=1)  stop("'order' must be a scalar.")
   if ( floor(order) != ceiling(order) )  stop("'order' must be an integer.")
   lwd=4; font=2; font.lab=2; las=1; cex.lab=2; cex.axis=1.8; cex.main=2
   subdivisions=1e4; tol=1e-13
   a0 = integrate( f, 0, 2*pi, subdivisions=subdivisions )$value / pi ;  a <- b <- NULL
   for ( i in 1:order ) {
      fcos = function(x) { match.fun( f )(x) * cos( i*x ) }
      fsin = function(x) { match.fun( f )(x) * sin( i*x ) }
      a = c( a, integrate( fcos, 0, 2*pi, subdivisions=subdivisions )$value/pi )
      b = c( b, integrate( fsin, 0, 2*pi, subdivisions=subdivisions )$value/pi )   }
   max.coef = max( abs(a0), abs(a), abs(b), tol*1e4 )
   if ( abs(a0) < max.coef*tol )   a0 = 0
   for ( i in 1:order ) {
      if ( abs(a[i]) < max.coef*tol )   a[i] = 0
      if ( abs(b[i]) < max.coef*tol )   b[i] = 0    }
   if (max.coef==0)   equation = "g(x) = 0"
   if (max.coef!=0)   {
      equation = "g(x) ="
      if ( a0 != 0 )  equation = paste( equation, as.character( prettyNum(a0/2,format="g",digits=3) ) )
      for ( i in 1:order )   {
         if ( a[i]>0 )  equation = paste( equation, " + ", 
                            as.character( prettyNum(a[i],format="g",digits=3) ), 
                            " cos ", ifelse(i>1.5,as.character(i),""), "x", sep="" )
         if ( a[i]<0 )  equation = paste( equation, " - ", 
                            as.character( prettyNum(-a[i],format="g",digits=3) ), 
                            " cos ", ifelse(i>1.5,as.character(i),""), "x", sep="" )
         if ( b[i]>0 )  equation = paste( equation, " + ", 
                            as.character( prettyNum(b[i],format="g",digits=3) ), 
                            " sin ", ifelse(i>1.5,as.character(i),""), "x", sep="" )
         if ( b[i]<0 )  equation = paste( equation, " - ",
                            as.character( prettyNum(-b[i],format="g",digits=3) ), 
                            " sin ", ifelse(i>1.5,as.character(i),""), "x", sep="" )   }      }  
   g = function(x){  output=a0/2 
       for ( i in 1:order )  output = output + a[i] * cos( i*x ) + b[i] * sin( i*x ) 
       return( output )   }
   x1 = -0.3*pi ;    x2 = 2.3*pi
   y1 = min( match.fun(f)( seq(0,2*pi,length.out=1000) ),
              match.fun(g)( seq(x1,x2,length.out=1000) ) )
   y2 = max( match.fun(f)( seq(0,2*pi,length.out=1000) ),
              match.fun(g)( seq(x1,x2,length.out=1000) ) )
   order.num = paste( order, "th", sep="" )
   if ( order %% 10 == 1 )  order.num = paste( order, "st", sep="" )
   if ( order %% 10 == 2 )  order.num = paste( order, "nd", sep="" )
   if ( order %% 10 == 3 )  order.num = paste( order, "rd", sep="" )
   plot( c(x1,x1,x2,x2), c(y1,y2,y1,y2), col="white", xlab="x", ylab="",
         main=paste( "f(x) = ", gsub(" ","",deparse(f)[3]), "    ", order.num, "order" ), xaxt="n",
         font=font, font.lab=font.lab, las=las, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, ... )
   abline( v=0, lty=3, lwd=3 ) ;    abline( v=2*pi, lty=3, lwd=3 )
   curve( f, 0, 2*pi, add=TRUE, xaxt="n", lwd=lwd*1.3, ... ) 
   curve( g, x1, x2, col="red", add=TRUE, xaxt="n", lwd=lwd, ... )    
   axis(1, at=c(0,pi/2,pi,3*pi/2,2*pi), 
         labels=c( 0, expression(paste(0.5,pi)), expression(pi), expression(paste(1.5,pi)), 
                      expression(paste(2,pi)) ), 
         font=font, font.lab=font.lab, cex.lab=cex.lab, cex.axis=cex.axis)
   return( list( constant=a0/2, cosine.coefficients=a, sine.coefficients=b, approximation=equation ) )
}
