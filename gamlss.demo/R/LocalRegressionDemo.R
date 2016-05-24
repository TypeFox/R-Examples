 demo.LocalRegression <- function(y=NULL, x=NULL,  span = 0.5,  position = trunc((n-1)/2) , 
                                 deg=1)  # w = rep(1, length(y)), 
 {
 #---------------------------------------------------
 # Simulate data
   if (is.null(y))
   {  
       n <- 100
      x <- seq(0, 1, length = n)*1.4
    set.seed(123)
    y <- 1.2 + .3*sin(5  * x) + rnorm(n) * 0.2
   }
      y <- y
      x <- if (is.null(x)) stop ("the x-variable must be set here") else x
   #}-------------------------------------------------------------
  span  <- span
  position <- trunc((n-1)/2)  
  deg <- 1
 if (interactive()) 
     {
     ps.smooth = function(panel)
      {
      span <- panel$span
      deg   <- panel$deg
      position <- panel$position
          LPOL(y, x, span=span, deg=deg, position=position) 
 
   panel
       }
         ps.panel = rp.control('Local Polinomial Smoother', size = c(400, 400),  span = 0.5)
       #  ps.panel = rp.control('Local Poly Smoother', size = c(400, 200),  span = 0.5)
          rp.slider(ps.panel, variable = span, from = 0.01, to = 2,  action = ps.smooth, resolution = 0.01, showvalue = TRUE, title = 'span')
           rp.doublebutton(ps.panel, variable = deg,  action = ps.smooth, initval = 1,  step = 1, range = c(0, 4), showvalue = TRUE, "degree")
         rp.slider(ps.panel, variable = position, from = 1, to = 100,  action = ps.smooth, resolution = 1,initval=49,  showvalue = TRUE, title = 'position')    
     }
 }
# #---------------------------------------------------------------
#--------------------------------------------------------------------

LPOL <- function(y, x, span = 0.5,  position = trunc((n-1)/2) , w = rep(1, length(y)), deg=1 )
{
                  n <- length(y)  
                 k <- trunc((n * span - 1)/2)  # the number obs below/above 
               o.x <- order(x)
             x.o <- x[o.x]
             y.o <- y[o.x]
            w.o <- w[o.x]
     if (any(w.o < 0)) 
        stop("weights should be positive")	
        lo <-  ifelse(position-k<1,1,position-k)
        up <- ifelse(position+k>n,n, position+k)
         lowx <-  x.o[lo] 
         upx   <-  x.o[up]
     plot(y.o~x.o,  type="n",  ylab = "y", xlab="x")        
     abline(v=lowx)
     abline(v=upx)
    polygon(c(lowx, upx, upx,lowx), c(min(y)-10,min(y)-100, max(y)+10, max(y)+10), col="grey90", lty=0)
   points(y.o[-position]~x.o[-position])
   points(y.o[position]~x.o[position], col="red", pch=18, cex=2)
  abline(v=x.o[position])
arrows(x.o[position], 0, x.o[position], 0.7, length =  0.25 ) 
yy <-  as.numeric(window(y.o, start = lo, end = up) )
xx <- as.numeric(window(x.o, start = lo, end = up))
   if (deg==0)
     { 
       fv <- rep(mean(yy), length(yy))	
       lines(fv~xx)	
       points(fv[(length(yy)-k)]~xx[(length(yy)-k)], col="green", cex=2, pch=4)
     }
     else
    {
      m1 <-lm(yy~poly(xx, deg))
      lines(fitted(m1)~xx)	
      points(fitted(m1)[(length(yy)-k)]~xx[(length(yy)-k)],  col="green", cex=2, pch=4)
    }
}
#--------------------------------------------------------------------
WLPOL<- function(y, x, sd = 0.5,  position = trunc((n-1)/2) , w = rep(1, length(y)), deg=1 )
{
                  n <- length(y)  
 #                k <- trunc((n * span - 1)/2)  # the number obs below/above 
               o.x <- order(x)
             x.o <- x[o.x]
             y.o <- y[o.x]
     #       w.o <- w[o.x]
     #if (any(w.o < 0)) 
      #  stop("weights should be positive")	
        
        #lo <-  ifelse(position-k<1,1,position-k)
        #up <- ifelse(position+k>n,n, position+k)
         #lowx <-  x.o[lo] 
         #upx   <-  x.o[up]
     plot(y.o~x.o,  type="n",  ylab = "y", xlab="x", ylim=c(0,1.71))        
   
    # xs <- seq(0,1.4,.05)
    nd <-  dnorm(x.o, mean=x[position], sd=.25)
    nds <- nd*.5
    lines(nds~x.o)
    # abline(v=lowx)
    # abline(v=upx)
    polygon(c(0,x.o, max(x)), c(0, nds, 0), col="grey90", lty=0)
   points(y.o[-position]~x.o[-position])
   points(y.o[position]~x.o[position], col="red", pch=18, cex=2)
  abline(v=x.o[position])
arrows(x.o[position], 0, x.o[position], 0.7, length =  0.25 ) 
#yy <-  as.numeric(window(y.o, start = lo, end = up) )
#xx <- as.numeric(window(x.o, start = lo, end = up))
   if (deg==0)
     { 
       fv <- rep(weighted.mean(y.o, w=nd), length(y.o))	
       lines(fv~x.o)	
       points(fv[1]~x[position], col="green", cex=2, pch=4)
     }
     else
    {
      m1 <-lm(y.o~poly(x.o, deg), weights=nd)
      lines(fitted(m1)~x.o)	
      fvp <- predict(m1, newdata=data.frame(x.o=x[position]))
      points(fvp~x[position],  col="green", cex=2, pch=4)
    }
}
#--------------------------------------------------------------------

  
  
   
