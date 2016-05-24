# $Id: fdhPlus.plot.R 156 2015-07-08 13:34:15Z b002961 $


dea.plot.fdhPlus <- function(x, y, param = 0.15, ...)  {

# Vaer sikker paa at x og y er matricer
x <- matrix(x)
y <- matrix(y)

# Saet parametrene low og high
if ( is.null(param) )  {
   delta <- .15
   low <- 1-delta
   high <- 1+delta
} else {
   if ( length(param) == 1 )  {
      low <- 1-param
      high <- 1+param
   } else {
      low <- param[1]
      high <- param[2]
   }
}



# Find FDH frontier
idx <- sort(x, index.return=TRUE)
effektive <- rep(NA, dim(x)[1])
j <- 1
prev <- idx$ix[1]
effektive[j] <- prev
for ( i in idx$ix )  {
   if ( y[i] > y[prev] )  {
   	j <- j+1
      effektive[j] <- i
      prev <- i
   }
}
# Frontier/rand og antal firms i randen
rand <- effektive[!is.na(effektive)]
x <- matrix(x[rand])
y <- matrix(y[rand])
# kun (x,y) paa FDH frontier er nu tilbage
sx <- sort(x,index.return=TRUE)
ix <- sx$ix
# crs-linjestykker starter i foerste soejles koordinat og 
# slutter i 2. soejles
x0 <- t(outer(c(low,high),x[ix]))
y0 <- t(outer(c(low,high),y[ix]))

if (FALSE)  {
points(x,y,pch=16,col="red")
text(x,y,1:dim(x)[1],adj=c(-.25,.75))
segments(x0[,1],y0[,1], x0[,2],y0[,2],col="lightgreen",lwd=1)
points(t(x0),t(y0),pch=c(2,0))
}

# text(x0[,1],y0[,1],1:dim(x0)[1],adj=c(-.75,.75))

# Det sidst tegnede linjestykke ender i (x_,y_)
segments(x0[1,1],0,  x0[1,1],y0[1,1] , ...)
x_ <- x0[1,1]
y_ <- y0[1,1]

nx <- dim(x0)[1]
for ( i in 1:(nx-1) )  {
   # Loeb alle de fdh effektive punkter igennem
   # Ved starten skal punkt vaere start eller punkt paa linjestykke 
   # der skal tegnes
# print(i)
# if ( i==15 ) break
   if ( x_ > x0[i,2] ) next; # vi er kommet laengere
   xk <- x_
   yk <- y_
   # Foerst ser vi om der er en linje over, dvs om der er en linje der
   # starter foer den nuvaerende slutter
   h <- 0
   while ( i+h+1<=nx && x0[i+h+1,1] < x0[i,2] )  {
      # Der er en linje foer det er slut, men er den ogsaa over?
      h <- h+1
      lambda <- (x0[i+h,1]-x0[i,1])/(x0[i,2]-x0[i,1])
      xk <- x0[i+h,1]
      yk <- (1-lambda)*y0[i,1] + lambda*y0[i,2]
      if ( y0[i+h,1] > yk ) break  # Linjen ligger over
   }
   if ( h > 0 && y0[i+h,1] > yk )  {
      # Der var en linje foer slut og den ligger oven over
      # Er vi allerede forbi det punkt?
      if ( x_ > xk || y_ > yk  )  next
      segments(x_,y_, xk,yk, ...) 
      segments(xk,yk,  x0[i+h,1],y0[i+h,1], ...)
      x_ <- x0[i+h,1]
      y_ <- y0[i+h,1]
      # vi er kommet til start paa ny linje
   }
      # Saa maa det vaere en lavere linje; find den forfra fordi den
      # maaske er sprunget over i soegen efter en hoejere
      
      # Den foerste linje maa vaere den lavere vi soeger
      # Er der en lavere linje?
   else if ( i+1<=nx && (x0[i+1,1] < x0[i,2] || y0[i+1,1] < y0[i,2]) )  {
         # der er en linje der starter lavere
         segments(x_,y_, x0[i,2],y0[i,2], ...)
         # Der er en lavere linje, find hvor den skal rammes vandret
         lambda <- (y0[i,2]-y0[i+1,1])/(y0[i+1,2]-y0[i+1,1])
         xk <- (1-lambda)*x0[i+1,1] + lambda*x0[i+1,2]
         yk <- y0[i,2]
         # Er der en hoejere linje foer xk?
         h <- 0
         while ( i+h<nx && x0[i+h+1,1] < xk )  {
            h <- h+1
            if ( y0[i+h,1] > yk ) break
         }
         if ( h > 0 && y0[i+h,1] > yk )  {
		      # der er en linje der starter over
            segments(x0[i,2],y0[i,2], x0[i+h,1], y0[i,2], ...)
            segments(x0[i+h,1], y0[i,2],  x0[i+h,1], y0[i+h,1], ...)
	         x_ <- x0[i+h,1]
	         y_ <- y0[i+h,1]
            # Vi er paa starten af en nyt linjesykke
	   } else {
            # Der maa saa vaere et stykke der starter under. 
            # Find det foerste der starter under
            lambda <- (y0[i,2]-y0[i+1,1])/(y0[i+1,2]-y0[i+1,1])
            xk <- (1-lambda)*x0[i+1,1] + lambda*x0[i+1,2]
            yk <- y0[i,2]
            h <- 1
            xkmin <- xk
            while ( i+h+1<=nx )  {
               # det er ikke til at vide hvor de sidste numre rammer 
               # fra neden
               h <- h+1
               lambda <- (y0[i,2]-y0[i+h,1])/(y0[i+h,2]-y0[i+h,1])
               if ( lambda < 0 || lambda > 1 ) break
               xk <- (1-lambda)*x0[i+h,1] + lambda*x0[i+h,2]
               if ( xk < xkmin )  {
                  xkmin <- xk
               }
            }
            segments(x0[i,2],y0[i,2], xkmin,yk, ...)
            x_ <- xkmin
            y_ <- yk
            # Vi er paa nyt linjestykke
         }
   } else {
      # Der var ingen linjestykker foer slut saa tegn linjen til slut
      segments(x_,y_, x0[i,2],y0[i,2], ...)
      x_ <- x0[i,2]
      y_ <- y0[i,2]
      # Find saa en vandret streg
      if ( i+1 > nx )  break
      if ( y0[i,2] > y0[i+1,1] )  {
         # Det naeste punkt starter lavere
         lambda <- (y0[i,2]-y0[i+1,1])/(y0[i+1,2]-y0[i+1,1])
         xk <- (1-lambda)*x0[i+1,1] + lambda*x0[i+1,2]
         yk <- y0[i,2]
         segments(x_,y_, xk,yk, ...)
         x_ <- xk
	   y_ <- yk
      } else {
	   # saa ma naeste punkt ligge over
         segments(x_,y_,  x0[i+1,1],y_, ...)
         segments(x0[i+1,1],y_,  x0[i+1,1],y0[i+1,1], ...)
         x_ <- x0[i+1,1]
         y_ <- y0[i+1,1]
      }
   }
}

segments(x_,y_, x0[nx,2],y0[nx,2], ...)
segments(x0[nx,2],y0[nx,2], 2*x0[nx,2],y0[nx,2], ...)


# }

}
