standard.ellipse <- function (x,y,confs=NULL,steps=5){


n <- length(x)

mx <- mean(x)
my <- mean(y)

# ------------------------------------------------------------------------------
## Original Parametric methods here
#varx <- var(x)
#vary <- var(y)
#
#sdx <- sd(x)
#sdy <- sd(y)
#
#r <- cov(x,y)/(sdx*sdy)
#
## Maths below taken from Batschelet. Circular Statistics in Biology
#A <- vary
#B <- -cov(x,y)
#C <- varx
#D <- (1-r^2)*A*C
#
#R <- ( (A-C)^2 + 4*B^2)^0.5
#a <- (2*D/(A+C-R))^0.5
#b <- (2*D/(A+C+R))^0.5
#theta <- atan(2*B/(A-C-R))
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Now do it using more succinct eigenvalue method
# ------------------------------------------------------------------------------

CM <- cov(cbind(x,y)) 
eig <- eigen(CM)

a <- sqrt(eig$values[1])
b <- sqrt(eig$values[2])

ac <- a*sqrt((n-1)/(n-2))
bc <- b*sqrt((n-1)/(n-2))

# NB this is the line causing odd ellipses to be drawn occasionally
#theta <- asin(eig$vectors[1,2]) # OLD LINE #
theta <- sign(CM[1,2]) * asin(abs(eig$vectors[2,1])) # NEW LINE 11/07/2011

SEA <- pi*a*b

psi <- seq(0,2*pi,steps*pi/180)

#xSEA <- numeric(steps+1)
#ySEA <- numeric(steps+1)

# ------------------------------------------------------------------------------
# calculate the coordinates of the standard ellipse
xtmp <- mx + a*cos(theta)*cos(psi) - b*sin(theta)*sin(psi)
ytmp <- my + a*sin(theta)*cos(psi) + b*cos(theta)*sin(psi)

tmp <- convexhull(xtmp,ytmp)

xSEA <- tmp$xcoords
ySEA <- tmp$ycoords

# ------------------------------------------------------------------------------
# now calculate the coordinates of the sample size corrected ellipse
xtmp <- mx + ac*cos(theta)*cos(psi) - bc*sin(theta)*sin(psi)
ytmp <- my + ac*sin(theta)*cos(psi) + bc*cos(theta)*sin(psi)

tmp <- convexhull(xtmp,ytmp)

xSEAc <- tmp$xcoords
ySEAc <- tmp$ycoords

# now calculate the confidence ellipse based on the provided levels

xCEA <- NULL
yCEA <- NULL
CEA <- NULL

if (!is.null(confs)){

  CEA <- numeric(length(confs))
  xCEA <- matrix(0,length(psi)+1,length(confs))
  yCEA <- xCEA
  v1 <- 2
  ct <- 0

  for (level in confs){

    ct <- ct + 1

    aCEA <- a * qf(level,v1,n-2)
    bCEA <- b * qf(level,v1,n-2)

    CEA[ct] <- pi * aCEA * bCEA

    xtmp <- mx + aCEA*cos(theta)*cos(psi) - bCEA*sin(theta)*sin(psi)
    ytmp <- my + aCEA*sin(theta)*cos(psi) + bCEA*cos(theta)*sin(psi)

    tmp <- convexhull(xtmp,ytmp)
    xCEA[,ct] <- tmp$xcoords
    yCEA[,ct] <- tmp$ycoords


  }
}

out <- list()
out$CEA <- CEA
out$SEA <- SEA
out$SEAc <- SEA * (n-1)/(n-2)
out$theta <- theta
out$confs <- confs
out$xCEA <- xCEA
out$yCEA <- yCEA
out$xSEA <- xSEA
out$ySEA <- ySEA
out$xSEAc <- xSEAc
out$ySEAc <- ySEAc
out$eccentricity <- sqrt(1-((b^2)/(a^2)))
out$a <- a
out$b <- b
#out$r <- r
out$ac <- ac
out$bc <- bc

return(out)
}