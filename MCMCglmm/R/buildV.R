"buildV"<-function(object, marginal=object$Random$formula, diag=TRUE, it=NULL, posterior="mean", ...){

  if(!is.null(marginal)){
    if(class(marginal)!="formula"){stop("marginal should be NULL or a formula")}
  }
  if(!is.null(it)){
    if(length(it)>1){stop("it should be an integer")}
    if(it>nrow(object$X) | it<1){stop("it should be less than or equal to the number of iterations")}
  }
  if(!is.null(posterior)){
    if(!posterior%in%c("distribution", "mean", "mode", "all")){
      stop("posterior argument must be either distribution, mean, mode or all")
    }
  }

  if(is.null(it)){
    if(posterior=="distribution"){
      it<-sample(1:nrow(object$X))
    }
    if(posterior=="all"){
      if(!diag){stop("diag must be equal to FALSE if posterior=='all'")}
    }
    if(posterior=="mean"){
      it<-1
      object$VCV[1,]<-colMeans(object$VCV)
    }
    if(posterior=="mode"){
      it<-1
      object$VCV[1,]<-posterior.mode(object$VCV, ...)
    }
  }

  rcomponents<-split.direct.sum(as.character(object$Random$formula)[2])

  if(length(rcomponents)!=length(object$Random$nrt)){stop("sorry - not implented for covu models")}

  mcomponents<-split.direct.sum(as.character(marginal)[2])

  if(any(mcomponents%in%rcomponents==FALSE)){stop("marginal formula does not correspond to model formula")}

  marginal<-rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt)

  if(diag){
    V<-matrix(0,nrow(object$Sol)*(posterior=="all")+(posterior!="all"),nrow(object$X))
  }else{
    V<-Diagonal(nrow(object$X),0)
  }

  cnt<-0
  cnt2<-0

  if(!is.null(object$Random$nfl)){
    for(j in 1:length(object$Random$nfl)){

     nfl<-object$Random$nfl[j]
     nrl<-object$Random$nrl[j]
     nat<-object$Random$nat[j]

     if(marginal[j]){

       if(posterior!="all"){
         Vtmp<-as(matrix(object$VCV[it,cnt2+1:(nfl^2)],nfl,nfl), "sparseMatrix")
       }else{
         Vtmp<-as(matrix(t(object$VCV[,cnt2+1:(nfl^2),drop=FALSE]),nrow=nfl), "sparseMatrix")
       }

       Ztmp<-object$Z[,cnt+1:(nrl*nfl),drop=FALSE]

       if(diag){
         if(posterior!="all"){
           if(nat==0){
             V<-V+apply(Ztmp,1,function(x){sum(Matrix(x,nrl,nfl)%*%Vtmp*x)})
           }else{
             if(any(colSums(Ztmp>0)>1)){
               V<-V+apply(Ztmp,1,function(x){sum(solve(object$ginverse[[nat]],Matrix(x,nrl,nfl))%*%Vtmp*x)})
             }else{
               DA<-Diagonal(x=diag(solve(object$ginverse[[nat]])))
               if(all(round(DA@x, 1e-10)==1)){
                 V<-V+apply(Ztmp,1,function(x){sum(Matrix(x,nrl,nfl)%*%Vtmp*x)})
               }else{
                 V<-V+apply(Ztmp,1,function(x){sum(DA%*%Matrix(x,nrl,nfl)%*%Vtmp*x)})
               }
             }
           }
         }else{
           if(nat==0){
             V<-V+apply(Ztmp,1,function(x){tapply(colSums(Matrix(x,nrl,nfl)%*%Vtmp*x), rep(1:nrow(object$Sol), each=nfl), sum)})
           }else{
             if(any(colSums(Ztmp>0)>1)){
               V<-V+apply(Ztmp,1,function(x){tapply(colSums((as(solve(object$ginverse[[nat]], Matrix(x,nrl,nfl)),"sparseMatrix")%*%Vtmp)*x), rep(1:nrow(object$Sol), each=nfl), sum)})
             }else{
               DA<-Diagonal(x=diag(solve(object$ginverse[[nat]])))
               if(all(round(DA@x, 1e-10)==1)){
                 V<-V+apply(Ztmp,1,function(x){tapply(colSums(Matrix(x,nrl,nfl)%*%Vtmp*x), rep(1:nrow(object$Sol), each=nfl), sum)})
               }else{
                 V<-V+apply(Ztmp,1,function(x){tapply(colSums((as(DA%*%Matrix(x,nrl,nfl),"sparseMatrix")%*%Vtmp)*x), rep(1:nrow(object$Sol), each=nfl), sum)})
               }
             }
           }
         }
       }else{
         if(nat==0){
           V<-V+Ztmp%*%kronecker(Vtmp, Diagonal(nrl))%*%t(Ztmp)
         }else{
           V<-V+Ztmp%*%kronecker(Vtmp, solve(object$ginverse[[nat]]))%*%t(Ztmp)
         }
       }
     }
     cnt<-cnt+(nrl*nfl)
     cnt2<-cnt2+nfl^2
   }
 }

  cnt<-0

  for(j in 1:length(object$Residual$nfl)){

    nfl<-object$Residual$nfl[j]
    nrl<-object$Residual$nrl[j]

    if(posterior!="all"){
      Vtmp<-matrix(object$VCV[it,cnt2+1:(nfl^2)],nfl,nfl)
    }else{
      Vtmp<-as(matrix(t(object$VCV[,cnt2+1:(nfl^2),drop=FALSE]),nrow=nfl), "sparseMatrix")
    }

    Ztmp<-object$ZR[,cnt+1:(nrl*nfl),drop=FALSE]

    if(posterior!="all"){
      if(diag){
        V<-V+apply(Ztmp,1,function(x){sum(Matrix(x,nrl,nfl)%*%Vtmp*x)})
      }else{
        V<-V+Ztmp%*%kronecker(Vtmp, Diagonal(nrl))%*%t(Ztmp)
      }
    }else{
      V<-V+apply(Ztmp,1,function(x){tapply(colSums(Matrix(x,nrl,nfl)%*%Vtmp*x), rep(1:nrow(object$Sol), each=nfl), sum)})
    }

    cnt<-cnt+(nrl*nfl)
    cnt2<-cnt2+nfl^2
  }
  
  return(V)
}
