#' Gibbs sampling of U and V
#' 
#' A Gibbs sampler for updating the multiplicative effect matrices U and V,
#' assuming they are the same across replicates. 
#' 
#' 
#' @usage rUV_rep_fc(E.T,U,V,rho,s2=1)
#' @param E.T Array of square residual relational matrix series with additive
#' effects and covariates subtracted out. The third dimension of the array is
#' for different replicates. Each slice of the array according to the third
#' dimension is a square residual relational matrix. 
#' @param U current value of U
#' @param V current value of V
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @author Peter Hoff, Yanjun He
#' @export rUV_rep_fc
rUV_rep_fc <-
  function(E.T,U,V,rho,s2=1)
  {
    Time<-dim(E.T)[3]
    R<-ncol(U)
    UV<-cbind(U,V)
    Suv<-solve(rwish(solve(diag(nrow=ncol(UV))+t(UV)%*%UV),2+nrow(UV)+ncol(UV)))
    
    Se<-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2<-mhalf(solve(Se))
    g<-iSe2[1,1] ; d<-iSe2[1,2]
    
    get.Er<-function(E,UVmr){
      return(E-UVmr)
    }
    get.Es<-function(Er,g,d){
      n<-sqrt(length(Er))
      return((g^2+d^2)*matrix(Er,n)+2*g*d*matrix(Er,n,byrow=T))
    }
    
    for(r in sample(1:R))
    {
      #Er.t<-Es.t<-array(dim=dim(E))
      #for (t in 1:Time){       
      #  Er.t[,,t]<- E[,,t] - U[,-r]%*%t(V[,-r]) ;  Es.t[,,t]<- (g^2+d^2)*Er.t[,,t]+2*g*d*t(Er.t[,,t])
      #}
      UVmr<-tcrossprod(U[,-r],V[,-r])
      Er.t<-apply(E.T,3,get.Er,UVmr)
      Es.t<-apply(Er.t,2,get.Es,g,d)
      
      vr<-V[,r]
      b0<- c(Suv[r,-r]%*%solve(Suv[-r,-r]))
      v0<- c(Suv[r,r] - b0%*%Suv[-r,r])
      m0<- cbind(U[,-r],V)%*%b0
      ssv<-max(sum(vr^2),1e-6)
      a<- Time*(g^2+d^2)*ssv+1/v0 ; c<- -2*Time*g*d/(a^2+a*2*Time*g*d* ssv)
      #Esv.vec<-apply(Es.t,1,sum)
      Esv.vec<-rowSums(Es.t)
      nEsv<-sqrt(length(Esv.vec))
      Esv<-matrix(Esv.vec,nEsv)%*%vr 
      m1<- Esv/a + c*vr*sum((Esv+m0/v0)*vr)  + m0/(a*v0)  
      #Vh<-sqrt(1/a)*diag(nrow(E))+(vr%*%t(vr))*(sqrt(1/a+ssv*c)-sqrt(1/a))/ssv 
      ah<-sqrt(1/a) ; bh<-(sqrt(1/a+ ssv*c)- sqrt(1/a) )/ssv ; e<-rnorm(nrow(E.T[,,1])) 
      U[,r]<- m1 + ah*e + bh*vr*sum(vr*e)
      

      
      ## update Vr 
      ur<-U[,r] 
      rv<-R+r
      b0<- c(Suv[rv,-rv]%*%solve(Suv[-rv,-rv]))
      v0<- c(Suv[rv,rv] - b0%*%Suv[-rv,rv])
      m0<- cbind(U,V[,-r])%*%b0
      ssu<-max(sum(ur^2),1e-6)
      a<- Time*(g^2+d^2)*ssu+1/v0 ; c<- -2*Time*g*d/(a^2+a*2*Time*g*d* ssu)
      tEsu<-matrix(Esv.vec,nEsv,byrow=T)%*%ur
      m1<-tEsu/a + c*ur*sum((tEsu+m0/v0)*ur)  + m0/(a*v0)  
      #Vh<-sqrt(1/a)*diag(nrow(E))+(ur%*%t(ur))*(sqrt(1/a+ssu*c)-sqrt(1/a))/ssu
      ah<-sqrt(1/a) ; bh<-(sqrt(1/a+ ssu*c)- sqrt(1/a) )/ssu ; e<-rnorm(nrow(E.T[,,1])) 
      V[,r]<- m1 + ah*e + bh*ur*sum(ur*e)
      ###
    }
     
    
    list(U=U,V=V)
  }
