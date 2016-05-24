minque <-
function(y, varcov, start=rep(0, length(k)), lower.bound=-Inf, restricted=TRUE)
{ 
	k=varcov; varRatio=start
  if(!isTRUE(restricted)) .NotYetImplemented()
  stopifnot(is.numeric(y) && is.list(k) && is.numeric(varRatio))
  nK=length(k)
  varRatio=c(rep(varRatio, length=nK),1)
  k[[nK+1L]]=diag(1,nrow(k[[1L]]))  
  
  LI=t(backsolve(chol(Reduce(`+`,mapply(`*`,k,varRatio,SIMPLIFY=FALSE))), k[[nK+1]]))
  LIk=lapply(lapply(k, `%*%`, LI), crossprod, LI)
  LIy=LI%*%y

  S=outer(LIk,LIk,function(a,b)sapply(mapply(`*`,a,b,SIMPLIFY=FALSE),sum))
  u=sapply(lapply(LIk,'%*%',LIy),crossprod,LIy)

  if(lower.bound<0 && is.infinite(lower.bound)){
    ans0=solve(S, u) ## FIXME: add error handler
  }else{
    # require(quadprog)
    qp.rslt=tryCatch(solve.QP(Dmat=crossprod(S), dvec=crossprod(S,u), Amat=rbind(diag(1,nK),0), bvec=rep(lower.bound,nK), meq=0L, factorized=FALSE), 
      error=function(e)  {
        if(e$message=='matrix D in quadratic function is not positive definite!'){
          solve.QP(Dmat=as.matrix(nearPD(crossprod(S))$mat), dvec=crossprod(S,u), Amat=rbind(diag(1,nK),0), bvec=rep(lower.bound,nK), meq=0L, factorized=FALSE) 
        }else{
          e
        }
      }
    )
    ans0=qp.rslt$solution
  }
  ans=ans0[-nK-1L]/ans0[nK+1L]
  ans
}
