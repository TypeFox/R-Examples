getsamples <-
function(mat,X,levelvec){
  n <- nrow(X)
  ind <- rep(TRUE,n)
  for (k in 1:nrow(mat)){
    var <- mat[k,"variable"]
    nfact <- length(levelvec[[var]])
    havenas <- which(is.na(X[,var]))
    if(length(havenas)>0) X[havenas,var] <- 1
    if(nfact==0){
      if(!is.inf((mat[k,"lower"]))) ind <- ind & (X[,var]> mat[k,"lower"])
      if(!is.inf((mat[k,"upper"]))) ind <- ind & (X[,var]<= mat[k,"upper"])
    }else{
      ind <- ind & (X[,var] %in% levelvec[[var]][which(dectobin(mat[k,"lower"],nl=nfact) >0.5)])
    }
    if(length(havenas)>0) ind[havenas] <- FALSE
  }
  if(any(is.na(ind))) ind[which(is.na(ind))] <- FALSE
  ind <- which(ind)
  return(ind)
}

