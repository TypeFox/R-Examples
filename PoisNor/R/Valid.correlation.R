Valid.correlation <-
function(no.pois, no.norm, lamvec){
  
  nPois = length(lamvec)
  
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


  if (no.pois!=nPois) { stop("Dimension of lamvec does not match the number of Poisson variables!\n") }
  if (sum(lamvec<0)>0){ stop("Lambda values cannnot be negative!\n") }
  
  if (no.norm<0) { stop("Number of normal variables cannnot be negative!\n") }
  if  (!is.wholenumber(no.norm)){ stop("Number of normal variables cannnot be fractional number!\n") }
  
  
  samples=100000
  u = runif(samples, 0, 1)
  X = rnorm(samples,0,1)
  Y = rnorm(samples,0,1)
  U = pnorm(X)
  
  maxmat = minmat = diag(NA,no.pois + no.norm)
  
  errorCount =0
  

  
  
  for (i in 1:(no.pois + no.norm) ){
    for (j in 1:(no.pois + no.norm)){
      
      if(j<=nPois){
      
        if (i!=j & i <=nPois ){ 
          
          maxcor=cor(qpois(u, lamvec[i]), qpois(u, lamvec[j]) )
          mincor=cor(qpois(u, lamvec[i]), qpois(1-u, lamvec[j]) )
          minmat [i,j] = mincor
          maxmat [i,j] = maxcor
          
        }  
        
        if (i > nPois) {
          
          Xpois = qpois(U,lamvec[j])
          
          max = cor(Xpois[order(Xpois)],Y[order(Y)])
          min = cor(Xpois[order(Xpois,decreasing=TRUE)],Y[order(Y)])
          minmat [i,j] = minmat [j,i] = min
          maxmat [i,j] = maxmat [j,i] =max        
        } 
        
      }else if(i> nPois & j > nPois){
        if (i!=j){
          minmat [i,j] = minmat [j,i] = -1
          maxmat [i,j] = maxmat [j,i] = 1  
        }
        
      }
      
      cat(".")
      
    }
    
  }
  cat("\n")

  return(list(min=minmat, max=maxmat))
}
