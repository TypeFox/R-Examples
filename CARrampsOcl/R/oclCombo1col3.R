oclCombo1col3 <-
function(a,b,b2,D, tausqy, tausqphi, By)
{
   if(!is.numeric(a) | !is.numeric(b) | !is.numeric(b2))
          stop("a, b, and b2 must be numeric matrices")

   na1 <- nrow(a)
   nb1 <- nrow(b)
   nb21 <- nrow(b2)
   nc1 <- length(tausqy)
   F1 <- ncol(tausqphi) + 1
   nab <- na1 * nb1 * nb21
   mresults <- rep(0, 2 * nab)

   out <- .C("oclCombo1col3", a=as.double(as.vector(t(a))), 
              b = as.double(as.vector(t(b))),
              b2 = as.double(as.vector(t(b2))),
              D = as.double(as.vector(t(D))),
              tausqy = as.double(tausqy) ,
              tausqphi = as.double( as.vector( t(tausqphi) )) ,
    By = as.double(By), results = as.double(mresults), 
               na1 = as.integer(na1), nb1 = as.integer(nb1),
               nb21 = as.integer(nb21), nc1 = as.integer(nc1), F1 = 
               as.integer(F1) )
   return(list( phimean = out$results[1:nab], phisd = sqrt(out$results[(nab+1):(2*nab)]/(nc1-1) )))
}

