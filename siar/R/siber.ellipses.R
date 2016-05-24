siber.ellipses <- function(x,y,group,R=10^4) {


# now loop through the data and calculate the ellipses
ngroups <- length(unique(group))



# split the isotope data based on group
spx <- split(x,group)
spy <- split(y,group)

SEA.B <- matrix(data=0,nrow=R,ncol=ngroups)

for (j in 1:ngroups){

  # fit a Multivariate Normal distribution to the data
  model <- bayesMVN(spx[[j]],spy[[j]],R = R)

  # alternatively you could fit two independent normal distributions
  #model <- bayestwoNorm(x,y,R=reps)

  Nobs <- nrow(model$b)

  for (i in 1:Nobs){

    estS <- model$S[i,]
    dim(estS) <- c(2,2)

    SEA.B[i,j] <- popSEA(estS)$SEA

  }

}


return(SEA.B)

}