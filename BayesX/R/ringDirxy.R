#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
## Project: BayesX
## Time-stamp: <[ringDirxy.R] by DSB Sam 28/02/2009 19:18 (GMT) on daniel@puc-home>
##
## Description:
## Private R code replacing maptools:::.ringDirxy.
##
## History:
## 28/02/2009   file creation
#####################################################################################

.ringDirxy <- function (xy)              # matrix of x and y coordinates as columns
{
    xy <- as.matrix(xy)
    nPoints <- nrow(xy)
    
    ## checks    
    stopifnot(is.numeric(xy),
              identical(ncol(xy), as.integer(2)),
              nPoints >= 3)
    
    ## index vector from the 2nd to the 2nd but last point
    inds <- seq_len(nPoints-2) + 1

    ## the start point and the middle and tail matrix
    start <- xy[1, ]                    # here the drop of dims is OK
    middle <- xy[inds, , drop=FALSE]    # here prevent drop of dims if there are only 3 points
    tail <- xy[inds + 1, , drop=FALSE]

    ## compute twice the signed areas
    areas <- (middle[, 1] - start[1]) * (tail[, 2] - start[2]) -
        (tail[, 1] - start[1]) * (middle[, 2] - start[2])
    ## and sum up to total twice signed area of the polygon
    total <- sum(areas)

    ## the sign then gives the direction
    if (total > 0)
    {
        return(as.integer(-1))          # counter-clockwise
    } else {
        return(as.integer(1))           # clockwise
    }
}


## closer analogue to original & C functions is here for documentation:

## .ringDirxy <- function (xy) 
## {
##     a <- xy[, 1]
##     b <- xy[, 2]
##     nvx <- length(b)

##     ## previous:
##     ## res <- .C("RFindCG", as.integer(nvx), as.double(a), as.double(b), 
##     ##     as.double(0), as.double(0), as.double(0), PACKAGE = "maptools")

##     ## now:
##     ## match args in RFindCG
##     n <- nvx
##     x <- a
##     y <- b
##     xc <- 0
##     yc <- 0
##     area <- 0

##     ## previous:
##     ## 	int i, nn;
##     ## 	tPointd *P;
##     ## 	tPointd CG;
##     ## 	double Areasum2;
##     ## 	nn = n[0];
##     ## 	P = (tPointd *) R_alloc(nn, sizeof(tPointd));
##     ## 	for (i=0; i<nn; i++) {
##     ## 		P[i][0] = x[i];
##     ## 		P[i][1] = y[i];
##     ## 	} 

##     ## now:
##     P <- cbind(x, y)

##     ## previous:
##     ## 	 FindCG(nn, P, CG, &Areasum2);
##     ## that is essentially
##     ## 	 for (i = 1; i < n-1; i++) {
##     ## 	        Centroid3( P[0], P[i], P[i+1], Cent3 );
##     ## 	        A2 =  Area2( P[0], P[i], P[i+1]);
##     ## 		CG[0] += A2 * Cent3[0];
##     ## 		CG[1] += A2 * Cent3[1];
##     ## 		Areasum2[0] += A2;
##     ## 	      }
##     ## 
##     ## now: (because we only need the Areasum2 here!)

##     ## Returns twice the signed area of the triangle determined by a,b,c,
##     ## positive if a,b,c are oriented ccw, and negative if cw.
##     twiceSignedArea <- function(a, b, c)
##     {
##         return((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]))
##     }

##     ## the loop from the 2nd to the (n-1)th point:
##     for(i in (seq_len(n-2) + 1))
##     {
##         area <- area + twiceSignedArea(a=P[1, ],
##                                        b=P[i, ],
##                                        c=P[i+1, ])
##     }
        
##     ## previous:
##     ## 	xc[0] = CG[0];
##     ## 	yc[0] = CG[1];
##     ## 	area[0] = Areasum2/2;
##     ## now: 
##     ## not necessary, as only the sign of area is relevant.
    
##     if (area > 0)
##     {
##         return(as.integer(-1))
##     } else {
##         return(as.integer(1))
##     }
## }
