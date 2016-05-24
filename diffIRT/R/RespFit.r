RespFit=function(object,order=2){
  if(object$conv!=0) stop("The solution in 'object' has not converged.")
  if(object$nr_fail!=0) cat("Warning: The solution provided in 'object' might be instable.\n")
  if(order==1) stop("The Mr test is not possible with only first-order moments, please use at least order=2")
  
    nit=object$nit
    N=object$N
    ai=object$par.log[1:nit]
    vi=object$par.log[(nit+1):(2*nit)]
    sd2A=exp(object$par.log[3*nit+1])
    sd2V=exp(object$par.log[3*nit+2])
    model=object$model
    nq=object$nq
    W=object$W
    A=object$A
    x=object$score
    constrain=object$constrain
      
    o=matrix(,1,1)                                                       # Vector for observed moments
    for(i in 1:order) o=cbind(o,fitObs(x,i))                         # Obsereved moments up to moment 'order'
    o=t(as.matrix(o[,-1]))
    pi.vec=fitPred(ai,vi,sd2A,sd2V,nit,model,A,W)                        # all expected moments = pi without dot in the matrix on the right upper page 1010; Heavy, in future try to omit its computation
    design=fitOrders(nit,nit)
    dim=fitCumChoose(nit,nit)
    r=fitCumChoose(nit,order)
    Tr=matrix(0,r,dim)
    Tr=cbind(0,matrix(.C("makeT",as.integer(nit),as.double(design),
            as.integer(dim),as.integer(r),as.double(Tr))[[5]],r,dim))    # matrix Tr from the upper right matrix and Eq 3 on p. 1010 of Maydeu - Olivares & Joe, 2005
    gamma=diag(c(pi.vec))-pi.vec%*%t(pi.vec)                             # gamma from the equation below Eq 2
    k=Tr%*%gamma%*%t(Tr)                                                 # ksi from Eq 3
    der=fitDeriv(ai,vi,sd2A,sd2V,nit,model,A,W)                          # calculate derivatives, expensive, maybe omit this calulation in the future
    if(model==2) der=der[,-((nit+1):(2*nit+2))]
    d=Tr%*%der                                                           # select the appropriate derivatives
    Mr=-999
    df=-999
    if(qr(d)$rank<ncol(d))
      cat("identification problems: matrix of first-order derivatives is not of full column rank\n")
    if(qr(d)$rank==ncol(d)) {
      dkd=t(d)%*%solve(k)%*%d
      sdkd=try(solve(dkd),silent=T)
      if(is.character(sdkd)) sdkd=try(solve(dkd+rnorm(ncol(dkd)^2,0,.0000001)),silent=T)
      if(!is.character(sdkd)){
        Cr=solve(k)-solve(k)%*%d%*%sdkd%*%t(d)%*%solve(k)
        Mr=N*t(t(o)-Tr%*%pi.vec)%*%Cr%*%(t(o)-Tr%*%pi.vec)
        info=fitCumChoose(nit,order)
        if(model==1) df=info-max(constrain[1:(2*nit)])-length(unique(constrain[(3*nit+1):(3*nit+2)]))-2
        if(model==2){
          x1=length(unique(constrain[1:nit]))
          x2=length(unique(constrain[(nit+1):(2*nit)]))
          df=info-max(x1,x2)-length(unique(constrain[(3*nit+1):(3*nit+2)]))-2
        }       
    } else cat("error: the weighted covariance matrix of the residuals is singular\n")
    }
    Z=rbind(t(Tr%*%pi.vec),(o))
    
    X=fitItems(nit,order)
    nms=c()
    for(i in 1:nrow(X)) nms[i]=paste(X[i,],collapse="")
    mat=t(Z)
    nrms=sqrt(N) * (mat[,2]- mat[,1])/sqrt((1/diag(Cr)))
    mat=cbind(round(mat,3),round(nrms,3))
    dimnames(mat)[[2]]=c('pred','obs','Z')
    mat=data.frame(mat,row.names=cbind(matrix(nms[-1])))
    respfitRES=list(model=model,Z=mat,Mr=Mr,df=df,order=order,nit=nit)
    class(respfitRES)<- "RespFit"
    return(respfitRES)
}
