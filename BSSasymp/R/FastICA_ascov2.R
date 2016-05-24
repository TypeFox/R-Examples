ASCOV_FastICAsym2 <- function(sdf, G, g, dg, supp=NULL, A=NULL, ...)
{
  p<-length(sdf)
  if(is.null(supp)) supp<-matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A<-diag(p)
  var_diag<-rep(0,p) 
 
  Ex4<-NULL
  Eg2<-NULL 
  Eg<-NULL
  Egx<-NULL
  Edg<-NULL  
  EG<-NULL  

  for(j in 1:p){
    Ex4[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j]<-(Ex4[j]-1)/4

    Eg2[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)^2}),
                      supp[j,1],supp[j,2],...)$value
    Eg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)}),
                     supp[j,1],supp[j,2],...)$value
    Egx[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)*x}),
                      supp[j,1],supp[j,2],...)$value
    Edg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*dg(x)}),
                      supp[j,1],supp[j,2],...)$value
    EG[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*G(x)}),
                      max(-20,supp[j,1]),min(20,supp[j,2]),...)$value
     
  }

   
  ASV<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j]<-(EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)+EG[j]^2*(Eg2[j]-Eg[j]^2-
2*Egx[j]*Edg[j]+Edg[j]^2))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
      }else ASV[i,j]<-var_diag[j]
    }
  }

  ASCOV<-diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i]<-(-EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)-EG[j]^2*(Eg2[j]-
Eg[j]^2-Egx[j]^2)+EG[i]*EG[j]*(Egx[i]*Egx[j]+Edg[i]*Edg[j]-Egx[i]*Edg[j]-
Egx[j]*Edg[i]))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
       
        ASCOV[(j-1)*p+i,(i-1)*p+j]<-ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  EMD<-sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W<-solve(A)
  W<-crossprod(diag(sign(rowMeans(W))),W)
  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)

}

ASCOV_FastICAsym_est2 <- function(X, G, g, dg, mixed=TRUE)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
 
  var_diag<-rep(0,p) 
  
  if(mixed){
    fI<-fICA(X,g,dg)
    W<-fI$W
  }else W<-diag(p)  

  X<-tcrossprod(sweep(X,2,colMeans(X)),W)   
  
  Ex4<-rep(0,p)
  Eg2<-rep(0,p) 
  Eg<-rep(0,p)
  Egx<-rep(0,p)
  Edg<-rep(0,p)  
  EG<-rep(0,p)  

  for(j in 1:p){
    Ex4[j]<-mean(X[,j]^4)
    var_diag[j]<-(Ex4[j]-1)/4

    Eg[j]<-mean(g(X[,j]))
    Eg2[j]<-mean(g(X[,j])^2)
    Egx[j]<-mean(g(X[,j])*X[,j])
    Edg[j]<-mean(dg(X[,j]))
    EG[j]<-mean(G(X[,j]))   
  }

  ASV<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j]<-(EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)+EG[j]^2*(Eg2[j]-Eg[j]^2-
2*Egx[j]*Edg[j]+Edg[j]^2))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
      }else ASV[i,j]<-var_diag[j]
    }
  }

  ASCOV<-diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i]<-(-EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)-EG[j]^2*(Eg2[j]-
Eg[j]^2-Egx[j]^2)+EG[i]*EG[j]*(Egx[i]*Egx[j]+Edg[i]*Edg[j]-Egx[i]*Edg[j]-
Egx[j]*Edg[i]))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
       
        ASCOV[(j-1)*p+i,(i-1)*p+j]<-ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A)
}



ASCOV_FastICAsym1 <- function(sdf, G, g, dg, supp=NULL, A=NULL, ...)
{
  p<-length(sdf)
  if(is.null(supp)) supp<-matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A<-diag(p)
  var_diag<-rep(0,p) 
 
  Ex4<-Eg2<-Eg<-Egx<-Edg<-EG<-sEGEgx_Edg<-NULL 

  EGn<-integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value
  for(j in 1:p){
    Ex4[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j]<-(Ex4[j]-1)/4

    Eg2[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)^2}),
                      supp[j,1],supp[j,2],...)$value
    Eg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)}),
                     supp[j,1],supp[j,2],...)$value
    Egx[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)*x}),
                      supp[j,1],supp[j,2],...)$value
    Edg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*dg(x)}),
                      supp[j,1],supp[j,2],...)$value
    EG[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*G(x)}),
                      max(-20,supp[j,1]),min(20,supp[j,2]),...)$value
    sEGEgx_Edg[j]<-sign(EG[j]*(Egx[j]-Edg[j])) 
  } 

  ASV<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j]<-(Eg2[i]-Eg[i]^2+Eg2[j]-Eg[j]^2-Egx[i]^2+Edg[j]*(Edg[j]-2*Egx[j]))/
               ((abs(Egx[i]-Edg[i])+abs(Egx[j]-Edg[j]))^2)
      }else ASV[i,j]<-var_diag[j]
    }
  }

  ASCOV<-diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i]<-(-(Eg2[i]-Eg[i]^2)-(Eg2[j]-Eg[j]^2)+Egx[i]^2+
       Egx[j]^2+sign(EG[i])*sign(EG[j])*(Egx[i]-Edg[i])*(Egx[j]-Edg[j]))/
       ((sign(EG[i])*(Egx[i]-Edg[i])+sign(EG[j])*(Egx[j]-Edg[j]))^2)
  
        ASCOV[(j-1)*p+i,(i-1)*p+j]<-ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  EMD<-sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W<-solve(A)
  W<-crossprod(diag(sign(rowMeans(W))),W)
  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)

}



ASCOV_FastICAsym21 <- function(sdf, G, g, dg, supp=NULL, A=NULL, ...)
{
  p<-length(sdf)
  if(is.null(supp)) supp<-matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A<-diag(p)
  var_diag<-rep(0,p) 
 
  Ex4<-NULL
  Eg2<-NULL 
  Eg<-NULL
  Egx<-NULL
  Edg<-NULL  
  EG<-NULL  

  EGn<-integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value
  for(j in 1:p){
    Ex4[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j]<-(Ex4[j]-1)/4

    Eg2[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)^2}),
                      supp[j,1],supp[j,2],...)$value
    Eg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)}),
                     supp[j,1],supp[j,2],...)$value
    Egx[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)*x}),
                      supp[j,1],supp[j,2],...)$value
    Edg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*dg(x)}),
                      supp[j,1],supp[j,2],...)$value
    EG[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*G(x)}),
                      max(-20,supp[j,1]),min(20,supp[j,2]),...)$value-EGn
     
  }

for(j in 1:p){
    Ex4[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j]<-(Ex4[j]-1)/4

    Eg2[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)^2}),
                      supp[j,1],supp[j,2])$value
    Eg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)}),
                     supp[j,1],supp[j,2])$value
    Egx[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*g(x)*x}),
                      supp[j,1],supp[j,2])$value
    Edg[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*dg(x)}),
                      supp[j,1],supp[j,2])$value
    EG[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*G(x)}),
                      max(-20,supp[j,1]),min(20,supp[j,2]))$value-EGn
     
  }
   

  ASV<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j]<-(EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)+EG[j]^2*(Eg2[j]-Eg[j]^2-
2*Egx[j]*Edg[j]+Edg[j]^2))/((abs(EG[i]*(Egx[i]-Edg[i]))+abs(EG[j]*(Egx[j]-Edg[j])))^2)
      }else ASV[i,j]<-var_diag[j]
    }
  }

  ASCOV<-diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i]<-(-EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)-EG[j]^2*(Eg2[j]-
Eg[j]^2-Egx[j]^2)+EG[i]*EG[j]*(Egx[i]-Edg[i])*(Egx[j]-Edg[j]))/
((abs(EG[i])*abs(Egx[i]-Edg[i])+abs(EG[j])*abs(Egx[j]-Edg[j]))^2)
       
        ASCOV[(j-1)*p+i,(i-1)*p+j]<-ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  EMD<-sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W<-solve(A)
  W<-crossprod(diag(sign(rowMeans(W))),W)
  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)

}



