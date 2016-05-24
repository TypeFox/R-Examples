"rIW"<-function(V,nu,fix=NULL, n=1, CM=NULL){
 
  if(is.matrix(V)==FALSE){stop("V must be a matrix")}

  # IMPORTANT #
  # Don't coerce scalar V's to a matrix - nu and V were swapped position from old version <1.10

  if(dim(V)[1]!=dim(V)[2]){stop("V must be square")}
  if(is.positive.definite(V, tol=1e-11)==FALSE){stop("V must be positive definite")}
  if(nu<=0){stop("nu must be greater than zero")}

  if(is.null(fix)){
    fix=-998
    CM<-1   
  }else{
    if(fix%%1!=0 | fix<1 | fix>dim(V)[1]){
      stop(paste("fix must be an integer between 1 and ", dim(V)[1]))
    }
    if(is.null(CM)){
      CM<-V[fix:dim(V)[1],fix:dim(V)[1],drop=FALSE]
    }else{
      if(is.matrix(CM)==FALSE){stop("CM must be a matrix")}
      if(dim(CM)[1]!=dim(CM)[2]){stop("CM must be square")}
      if(is.positive.definite(CM, tol=1e-11)==FALSE){stop("CM must be positive definite")}
      if(dim(CM)[1]!=(dim(V)[1]-fix+1)){stop(paste("expecting a ", (dim(V)[1]-fix+1), "x", (dim(V)[1]-fix+1), " matrix for CM", sep=""))}		   
    }
    if(sum(CM!=0)>dim(CM)[1]){stop("sorry - matrices to be conditioned on must be diagonal")}
  }

  if(fix==1){
    if(n==1){
      return(V)
    }else{
      return(t(matrix(rep(V,n), dim(V)[1]^2, n)))
    }
  }else{
    output<-.C("rIW",
      as.double(nu),
      as.double(c(solve(V*nu))),
      as.double(c(CM)),			   
      as.integer(dim(V)[1]),
      as.integer(fix-1),
      as.integer(n),
      as.double(matrix(0,n,dim(V)[1]^2))
    )
    if(n==1){
      return(matrix(output[[7]], dim(V)[1], dim(V)[1]))
    }else{
      return(t(matrix(output[[7]], dim(V)[1]^2, n)))
    }
  }
}
