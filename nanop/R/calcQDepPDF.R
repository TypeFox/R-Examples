calcQDepPDF <- function(nanop=NA, dr=.1, minR=1, maxR=20, dQ=.01, minQ=1, maxQ=20,
                        verbose=0, subdivisions = 100, order=1000,
                        rel.tol=.Machine$double.eps^.7,
                        addNoise=FALSE,
                        noiseFun=NA, 
						totalScattParams=list(), 
                        preTotalScat=NA, ...) {
 
  r <- seq(minR, maxR, by = dr)
  res <- vector(length=length(r))
  tnanop <- t(nanop)
  
  sigma=totalScattParams$sigma
  n=totalScattParams$n
  delta=totalScattParams$delta
  kind=totalScattParams$kind 
  dr = totalScattParams$dr
  del = totalScattParams$delta
  eps=totalScattParams$eps
  type=totalScattParams$type
  scatterLength=totalScattParams$scatterLength
  scatterFactor=totalScattParams$scatterFactor	
  
  if(is.na(preTotalScat[[1]][1])) 			   
    totalScatt <- calcTotalScatt(nanop, dQ=dQ, minQ=minQ, maxQ=maxQ,
                                 scatterFactor=scatterFactor, scatterLength=scatterLength,  
                                 sigma=sigma, n=n, delta=delta, kind=kind, type=type, 
						         dr = dr,  del = del, eps=eps)
  else
    totalScatt <- preTotalScat 
      
    
  for(i in 1:length(r)) {   
    res[i] <- distrEx::GLIntegrate(f=calcQDepPDFAux, 
	                      lower=minQ, upper=maxQ,
						  order=order, rel.tol=rel.tol,
						  r=r[i],  
						  tnanop=as.vector(tnanop),
						  minR=minR,
						  subdivisions=subdivisions, 
                          totalScatt = totalScatt,
                          addNoise=addNoise, 
                          noiseFun=noiseFun, ...) * (2/pi)
    if(verbose > 0 && ((i %% verbose) == 0))
      cat("Finished computing r=", r[i], "\n")
  }
  
  list(r=r, gr=res)
}

calcQDepPDFAux <- function(Q, r, 
                           totalScatt = NA,
                           addNoise=FALSE,
                           noiseFun=NA, ...) {

  fn1  <- function(res) noiseFun(res, ...)
  totalScattVec <- approx(totalScatt$Q, totalScatt$gQ, Q)$y

  if(addNoise)
    totalScattVec <- noiseFun(totalScattVec, ...)

  res <- .C("calcQDepPDF",
              res = as.double(Q),
              Q = as.double(Q),
              r = as.double(r), 
              len = as.integer(length(Q)),
              totalScattVec = as.double(totalScattVec), 
              PACKAGE="nanop")$res

  res
}
