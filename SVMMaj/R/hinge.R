getHinge <- function(hinge = 'quadratic', hingek = 3, eps=1e-8) {
#=========================================================
#GENERATE HINGE ERROR FUNCTION

#DETERMINE THE HINGE TERM AND MAJOR PARAMETERS
    hingechoice <- switch(hinge,
      absolute = function(q,y)        {
          z    <- q * y
          m    <- abs( 1 - z )
          m    <- m * ( m > eps ) + eps * ( m <= eps )
          a    <- .25 * m ^ -1
          b    <- y * ( a + .25 )
          c    <- a + .5 + .25 * m
          loss <- (1 > z)*( 1 - z )
          list( a=a , b=b , c=c , loss=loss ,call=match.call() )     },

      quadratic = function(q,y)     {
          z    <- q * y
          m    <- z * ( z > 1 ) + ( 1 >= z )
          a    <- 1
          b    <- y * m
          c    <- m ^ 2
          loss <- (1 > z)*( 1 - z )^2
          list( a=a , b=b , c=c , loss=loss ,call=match.call() )     },
  
      huber = function(q,y)        {
          z    <- q * y
          m    <- ( z < -hingek ) * (z + hingek)
          o    <- ( z >  1 ) * (z - 1)
          a    <- .5 * (hingek + 1)^-1
          b    <- y * a * ( 1 + m + o )
          c    <- y * b * ( 1 + m + o ) - m
          loss <- a*q^2 - 2*b*q + c
          list( a=a , b=b , c=c , loss=loss ,call=match.call() )     },

      logit = function(q,y)     {
          z    <- q * y
          m    <- ( 1 + exp(z) ) ^-1
          l    <- log( 1 + exp(-z)  )
              l.inf    <- is.infinite(l)
              l[l.inf] = log(1+exp(z[l.inf]))-z[l.inf]
          a    <- .125
          b    <- y * a * ( 4 * m + z )
          c    <- z * a * ( 8 * m + z )  + l
          list( a=a , b=b , c=c , loss=l ,call=match.call() )   },
          
#      logit2 = function(q,y)    {
#          z    <- q * y
#          m    <- ( 1 + exp(z) ) ^-1
#          l    <- log( 1 + exp(-z) )
#              l.inf    <- is.infinite(l)
#              l[l.inf] = log(1+exp(z[l.inf]))-z[l.inf]
#          z.a  <- pmax(abs(z), eps)
#          a    <- .25*(1-exp(-z.a))/( z.a*   (  1+exp(-z.a)))
#          b    <- y   * ( a * z + .5 * m  )
#          c    <- z * ( 2 * y * b - a * z ) + l
#          list( a=a , b=b, c=c , loss=l , call=match.call() )  },

#       probit = function(q,y) {
#        	z    <- q * y
#        	M    <- pnorm(z)
#        	l    <- -log(M)
#        	mr   <- dnorm(z)/M
#        	mr[M==0] <- -z[M==0]
#        	l[M==0]  <- .5*z[M==0]^2
#        	a    <- .5
#        	b    <- y*(z + .5*mr)
#        	c    <- z*(2*b*y - z) + l
#        	list( a=a, b=b, c=c, loss=l, call=match.call())   },

#       regression =	function(q,y){
#          qa   <- q - y
#          m    <- pmax( abs( abs(qa) - hingek ) , eps )
#          z    <- abs( qa ) - 2 * hingek
#		      a    <- .5 / ( (z > 0) * abs(qa)  +  (z <= 0) * 2 * m )
#		      b    <- (z < 0) * a * abs(hingek - m) * sign(qa)
#		      c    <- (z < 0) * a * ( hingek - m )^2 + (z > 0) * ( .5 * abs(qa) - hingek )
#		      loss <- pmax( abs( qa ) - hingek, 0 )
#		      list( a=a , b=b , c=c , loss=loss , call=match.call() ) }
    )
  class(hingechoice)          <- 'hinge'
  attr(hingechoice,'hinge')   <- hinge
  attr(hingechoice,'fixed.a') <- (hinge!='absolute' && hinge!='regression')
  if(hinge=='huber' || hinge=='regression') attr(hingechoice,'hingek') <-hingek
  if(hinge=='absolute' || hinge=='regression') attr(hingechoice,'eps') <- eps
  return(hingechoice)
}

#====================================================================
#====================================================================
#====================================================================

print.hinge <- function(x,...) {
	cat('Hinge function\n','  Hinge type:',attr(x,'hinge'),'\n')
	if(!is.null(attr(x,'hingek')) )
	  cat('  Hinge parameters:\n','k = ',attr(x,'hingek'),'\n')
}


#====================================================================
#====================================================================
#====================================================================

plot.hinge <- function(x,y=1,z=NULL,...) {
	object<-x
	#DEFINE PLOT RANGE
	if(!is.null(z))  scale <- max(abs(z)*1.2,5)
	else             scale <- 5
	#CALCULATE LOSS FUNCTION
	xs    <- seq(-1,1,by=.02) * scale
	loss <- object(xs,y)$loss
	plot(xs,loss,'l',main='Plot of hinge error',xlab='q',ylab='hinge value',...)
	#CALCULATE, IF SPECIFIED MAJORIZATION FUNCTION
	if(!is.null(z)){
	  hingechoice <- object(z,y)
	  hingeVal    <- hingechoice$a * xs^2 - 2 * hingechoice$b * xs + hingechoice$c
	  lines(xs,hingeVal,col='blue')
	  points(z,hingechoice$loss)
	}
}

