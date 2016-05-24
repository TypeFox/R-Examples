make.del<-function(pars)
  {
    # Takes a parameter vector with log of diagonal elements first
    # and then elements above diagonal in column-descending order
    # and fills in the zeros to make an upper-triangular matrix
    k<-floor((-1+sqrt(1+8*length(pars)))/2)
    mymatrix<-diag(exp(pars[1:k]))
    pars<-pars[-(1:k)]
    if (k>1){
      for(i in 2:k){
        mymatrix[1:(i-1),i]<-pars[1:(i-1)]
        pars<-pars[-(1:(i-1))]
      }
    }
    mymatrix
  }

