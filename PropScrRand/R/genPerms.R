genPerms <-
function(n, n1, nPerms){

  if(n < 1)            stop('n must be a positive integer')
  if(n1 < 1 | n1 >= n) stop('n1 is number to be treated: need 0 < n1 < n')

  ## equal allocation
  if(n1 == n/2){
  
    ## reset nPerms if choose(n, n1) is small
    if(nPerms >= choose(n, n1)/2){
      nPerms = choose(n, n1)
      warning(paste('nPerms too large, changing to choose(',
                    n,', ',n1,')',sep=''))
    }
    
    if(nPerms == choose(n, n1)){
      treatIndex = combn(n, n1)
    }else{
      tmat1 = replicate(n=nPerms, 
                        expr=sample(rep(0:1, length.out=n), replace=FALSE))
      tmat2 = unique(cbind(tmat1, 1-tmat1), MARGIN=2)

      while(ncol(tmat2) < 2*nPerms){
        needed = nPerms - ncol(tmat2)/2
        tmat3 = replicate(n=needed, 
                          expr=sample(rep(0:1, length.out=n), replace=FALSE))
        tmat4 = cbind(tmat2, tmat3, 1-tmat3)
        tmat2 = unique(tmat4, MARGIN=2)
      }
      
      ## get treatment indices
      treatIndex = sapply(1:ncol(tmat2), function(i){ which(tmat2[,i] == 1) })
    }
  }
  
  ## unequal allocation
  else{
    
    ## reset nPerms if choose(n, n1) is small
    if(nPerms >= choose(n, n1)){
      nPerms = choose(n, n1)
      warning(paste('nPerms too large, changing to choose(',
                    n,', ',n1,')',sep=''))
    }
  
    if(nPerms == choose(n, n1)){
      treatIndex = combn(n, n1)
    }else{
      tmat1 = replicate(n=nPerms, 
                        expr=sample(c(rep(1,n1),rep(0,n-n1)), replace=FALSE))
      tmat2 = unique(tmat1, MARGIN=2)
      while(ncol(tmat2) < nPerms){
        tmat3 = replicate(n=nPerms, 
                          expr=sample(c(rep(1,n1),rep(0,n-n1)),replace=FALSE))
        tmat4 = cbind(tmat2,tmat3)
        tmat2 = unique(tmat4, MARGIN=2)
      }
      tmat = tmat2[,1:nPerms]

      treatIndex = sapply(1:ncol(tmat), function(i){ which(tmat[,i] == 1) })
    }
  }
  
  return(treatIndex)
}
