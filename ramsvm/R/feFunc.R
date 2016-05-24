feFunc <- function(varTemp, 
                   nobs, 
                   k, 
                   fold, 
                   len, 
                   labidx, 
                   sam,
                   y.train,
                   xinner,
                   weight,
                   warm,
                   yyi,
                   ytrain,
                   templambda,
                   my,
                   gamma){

  kdouble <- k

  warm3 <- matrix(data = 0.0, nrow = nobs, ncol = k)
     
  index <- numeric(0)

  for( i in 1L:k ) {

    if( varTemp < fold ) {
      rng <- ((len[i]*(varTemp-1)+1)) : (len[i]*varTemp)
    } else {
      rng <- ((len[i]*(varTemp-1)+1)) : length(labidx[i])
    }
    select.index <- labidx[[i]][(sam[[i]])[rng]]  

    index <- c(index,select.index)
  }

  y.train.temp <- y.train[index]
  nobs.temp <- length(y.train.temp)
  xinner.temp <- xinner[index,index]
  weight.temp <- weight[index]
  alpha_ij.temp <- warm[index,]
  yyi.temp <- yyi[index,]
  nobsdouble.temp <- as.double(nobs.temp)
  ytrain.temp <- ytrain[index]

  alpha_yi.temp <- rep(0,nobs.temp)

  for( zz in 1L:nobs.temp ) {
    alpha_yi.temp[zz] <- alpha_ij.temp[zz,y.train.temp[zz]]
  }

  erci.temp <- -diag(xinner.temp)/2/nobs.temp/templambda

  aa <- .C("local_alpha_update", 
           as.vector(alpha_ij.temp),  
           as.vector(alpha_yi.temp),
           as.vector(my),  
           as.vector(yyi.temp),  
           as.vector(xinner.temp),
           as.double(templambda),  
           as.vector(weight.temp),  
           as.integer(nobs.temp),
           as.double(nobsdouble.temp),  
           as.integer(k),  
           as.double(kdouble),
           as.vector(erci.temp),  
           as.double(gamma),  
           as.vector(ytrain.temp), 
           outalpha_ij = as.vector(rep(0,nobs.temp*k)),
           Package = "ramsvm")

  warm3[index,] <- matrix(aa$outalpha_ij,nobs.temp,k)

  return(warm3)

}

