#
#  Copyright (C) 1997-2010 Friedrich Leisch
#  $Id: mlbench-regression.R 4612 2010-10-08 09:51:20Z leisch $
#

mlbench.friedman1 <- function(n, sd=1){

  x <- matrix(runif(10*n),ncol=10)

  y <- 10 * sin(pi * x[,1] * x[,2])
  y <- y + 20 * ( x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]

  if(sd>0){
    y <- y + rnorm(n, sd=sd)
  }

  list(x=x, y=y)
}

mlbench.friedman2 <- function(n, sd=125){

  x <- cbind(runif(n,min=0,max=100),
	     runif(n,min=40*pi,max=560*pi),
	     runif(n,min=0,max=1),
	     runif(n,min=1,max=11))

  y <- sqrt(x[,1]^2 + (x[,2]*x[,3] - 1/(x[,2]*x[,4]))^2)

  if(sd>0){
    y <- y + rnorm(n, sd=sd)
  }

  list(x=x, y=y)
}

mlbench.friedman3 <- function(n, sd=0.1){

  x <- cbind(runif(n,min=0,max=100),
	     runif(n,min=40*pi,max=560*pi),
	     runif(n,min=0,max=1),
	     runif(n,min=1,max=11))

  y <- atan( (x[,2]*x[,3] - 1/(x[,2]*x[,4])) / x[,1] )

  if(sd>0){
    y <- y + rnorm(n, sd=sd)
  }

  list(x=x, y=y)
}

mlbench.peak <- function(n, d=20)
  {
    metro <- numeric(n)
    y <- numeric(n)
    x <- matrix(0, nrow=n, ncol=d)
    for (ndata in 1:n)
      {
        radius <- runif(1, min=0, max=3)
        x[ndata,] <- rnorm(d)
        metro[ndata] <- sqrt(sum(x[ndata,]^2))
        x[ndata,] <- radius * (x[ndata,]/metro[ndata])
        y[ndata] <- 25 * exp(-0.5* radius^2)
      }
    list(x=x, y=y)
  }






          
        




