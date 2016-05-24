####################################
#This functions obtains the prior  #
#distribution speficied by the user#
####################################

which.dist <- function(dist)
  {
  split1 <- unlist(strsplit(dist,  "\\(" ))
  name.dist <- paste("d",split1[1],sep="")
  par.dist <- unlist(strsplit(split1,  "\\)" ))[2]

  a <- as.numeric(unlist(strsplit(par.dist,  "," ))[1])
  b <- as.numeric(unlist(strsplit(par.dist,  "," ))[2])
    
  ind <- switch(name.dist, dunif = , i=1, dgamma =, i = 2, dexp = , i = 3,
         dnorm = , i = 4, dt = , i=5, dweibull = ,i=6, df = ,i=7, dchisq = , i=8,
         dcauchy = , i=9, dlnorm = , i=10)

  if(is.na(b)){b=0}

  return(c(ind, a, b))
  }
   
    

