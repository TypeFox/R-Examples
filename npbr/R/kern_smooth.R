kern_smooth<-function(xtab, ytab, x, h, method="u", technique="noh")
{ 
  stopifnot(method%in%c("u","m","mc")) 
  stopifnot(technique%in%c("pr","noh")) 
    
  #internal parameters
  consX<-unique(sort(c(xtab,seq(min(xtab),max(xtab),length=201)))) # to ensure the obedience to the shape constraints
  Cstable=10000 # stability parameter
  
  n<-length(xtab)
  ndata<-length(xtab) # number of data points
  ncons<-length(consX) # number of the points where the constraints are imposed
  
  # sorting step to use the Priestly-Chao estimator
  oind<-order(xtab)
  xtab<-xtab[oind]
  ytab<-ytab[oind]
  
  if (technique=="noh")
  { 
  # sorting step to use the Priestly-Chao estimator
  oind<-order(xtab)
  xtab<-xtab[oind]
  ytab<-ytab[oind]
  
  xtab2<-xtab[-1]
  ytab2<-ytab[-1]
  DX<-diff(xtab)
  
  
  r_end <- max(xtab)
  l_end <- min(xtab)
  
  opt_coef<-(pnorm((r_end-xtab2)/h)-pnorm((l_end-xtab2)/h))*DX*ytab2
  
  C0<-.fitMat(xtab, ytab, xtab, h, type=0) # only envelop the data
  C1<-.fitMat(xtab, ytab, consX, h, type=1) # monotonicity constraints at the given points (including the data points)
  C2<-.fitMat(xtab, ytab, consX, h, type=2) # convexity constraints at the given points (including the data points)
  
  obj<-c(opt_coef,-opt_coef)
  
  if(method=="u")
  {
   mat_temp<-C0
   mat<-rbind(cbind(mat_temp,-mat_temp),rep(1,length(obj)))
   rhs<-c(ytab,Cstable)
   dir <- c(rep(">=",ndata),"<=")
  }
  
  
  if (method=="m")
  {
    mat_temp<-rbind(C0,C1)
    mat<-rbind(cbind(mat_temp,-mat_temp),rep(1,length(obj)))
    rhs<-c(ytab,rep(0,ncons),Cstable)
    dir <- c(rep(">=",ndata+ncons),"<=")
  }
  
  if (method=="mc")
  {       
    mat_temp<-rbind(C0,C1,C2)
    mat<-rbind(cbind(mat_temp,-mat_temp),rep(1,length(obj)))
    rhs<-c(ytab,rep(0,2*ncons),Cstable)
    dir <- c(rep(">=",ndata+ncons),rep("<=",ncons),"<=")
  }
    
  Sol<-Rglpk_solve_LP(obj, mat, dir, rhs, types=NULL, max=FALSE)
  OPT_temp<-Sol$sol
  # p+ and p-
  OPT<-OPT_temp[1:length(opt_coef)]-OPT_temp[(length(opt_coef)+1):(length(opt_coef)*2)]
  
  fitt<-.fitMat(xtab,ytab,x,h,0) %*% OPT
  }
  
  if (technique=="pr")
  {
    A <- .fitMat_nw(xtab, ytab, xtab, h, type=0) # only envelop the data
    A.deriv.1 <- .fitMat_nw(xtab, ytab, consX, h, type=1) # monotonicity constraints at the given points (including the data points)
    A.deriv.2 <- .fitMat_nw(xtab, ytab, consX, h, type=2) # convexity constraints at the given points (including the data points)
    
    if (method=="u")
    {
      p.hat <- solve.QP(Dmat=diag(n),
                        dvec=rep(1,n),
                        Amat=t(A),
                        bvec=ytab,
                        meq=0)$solution
    }
    
    if (method=="m")
    {
      p.hat <- solve.QP(Dmat=diag(n),
                        dvec=rep(1,n),
                        Amat=cbind(t(A),t(A.deriv.1)),
                        bvec=c(ytab,rep(0,ncons)),
                        meq=0)$solution
    }                  
    
    if (method=="mc")
    {
      p.hat <- solve.QP(Dmat=diag(n),
                        dvec=rep(1,n),
                        Amat=cbind(t(A),t(A.deriv.1),-t(A.deriv.2)),
                        bvec=c(ytab,rep(0,ncons),rep(0,ncons)),
                        meq=0)$solution
    }               
    
    fitt <- .fitMat_nw(xtab, ytab, x, h, type=0)%*%p.hat
   } 
  
  return(as.vector(fitt))
  
}

.dderi<-function(x){return(-x*dnorm(x))}

.d2deri<-function(x){return((x^2-1)*dnorm(x))}

.fitMat<-function(X, Y, consX, h, type=0)
{
  # This function constructs design matrix for estimation based on the Priestly-Chao estimator
  # X : vector of x-points of the data
  # Y : vector of y-points of the data
  # consX : vector of X-points on which the shape constraints are imposed 
  # h : bandwidth
  
  # type=0 for the function itself
  # type=1 for its first derivatives
  # type=2 for its second derivatives
  
  n<-length(X)
  stopifnot(sum((order(X)-1:n)^2)==0) # the data should be sorted so that X_1 \leq ... \leq X_n
  
  X2<-X[-1]
  Y2<-Y[-1]
  
  DX<-diff(X)
  A<-dnorm(outer(consX,X2,'-')/h)/h
  B<-rep(1,length(consX)) %*% t(Y2*DX)
  
  if (type==0) {return(A*B)}
  
  if (type==1) 
  {
    A2<-.dderi(outer(consX,X2,'-')/h)/h^2
    return(A2*B)
  }
  
  if (type==2)
  {
    A3<-.d2deri(outer(consX,X2,'-')/h)/h^3
    return(A3*B)
  }
}


.fitMat_nw<-function(X, Y, consX, h, type=0)
{
  # This function constructs design matrix for estimation based on the Nadaraya-Watson estimator
  # X : vector of x-points of the data
  # Y : vector of y-points of the data
  # consX : vector of X-points on which the shape constraints are imposed 
  # h : bandwidth
  
  # type=0 for the function itself
  # type=1 for its first derivatives
  # type=2 for its second derivatives
  
  
  A<-dnorm(outer(consX,X,'-')/h)/h
  ksum<-apply(A,1,sum)
  B<-(1/ksum)*rep(1,length(consX)) %*% t(Y)
  
  A2<-.dderi(outer(consX,X,'-')/h)/h^2
  kdsum<-apply(A2,1,sum)
  
  B2_1<-(1/ksum)*rep(1,length(consX)) %*% t(Y)
  B2_2<-(kdsum/ksum^2)*rep(1,length(consX)) %*% t(Y)
  
  A3<-.d2deri(outer(consX,X,'-')/h)/h^3
  kd2sum<-apply(A3,1,sum)
  
  B3_1<-(kd2sum/ksum^2)*rep(1,length(consX)) %*% t(Y) 
  B3_2<-(kdsum^2/ksum^3)*rep(1,length(consX)) %*% t(Y) 
  
  if (type==0) {return(A*B)}
  
  if (type==1) 
  {
    return(A2*B2_1-A*B2_2)
  }
  
  if (type==2)
  {
    alpha<-A3*B2_1-A2*B2_2
    beta<-A2*B2_2+A*B3_1-2*A*B3_2
    return(alpha-beta)
  }
}

