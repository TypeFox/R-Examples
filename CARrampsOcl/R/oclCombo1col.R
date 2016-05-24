oclCombo1col <- function(a,b,D, tausqy, tausqphi, By)
{
   if(!is.numeric(a) | !is.numeric(b)  )
          stop("a and b must be numeric matrices")

   na1 <- nrow(a)
   nb1 <- nrow(b)
   nc1 <- length(tausqy)
   F1 <- ncol(tausqphi) + 1
   nab <- na1 * nb1
   results <- rep(0, 2 * nab)

   out <- .C("oclCombo1col", a=as.double(as.vector(t(a))), 
              b = as.double(as.vector(t(b))),
              D = as.double(as.vector(t(D))),
              tausqy = as.double(tausqy) ,
              tausqphi = as.double( as.vector( t(tausqphi) )) ,
    By = as.double(By), results = as.double(results),
               na1 = as.integer(na1),
               nb1 = as.integer(nb1), nc1 = as.integer(nc1), F1 = 
               as.integer(F1) )
   return(list( phimean = out$results[1:nab], phisd = sqrt(out$results[(nab+1):(2*nab)]/(nc1-1) )))
}

