#####functions to generate the data
##### 5-22-07

######################################## orthornoral basis functions
##orthonormal basis  {Phi_1,..., Phi_{M}}  on [0,1];

######
## (i) fourier basis on [0,1]: sqrt(2)sin(k*pi*t), k=1,2,..,: "sin"
Basis.Sin<-function(t,M){
##para: t--evaluation point; M:dimension
##return: basis functions evaluated at t; M by 1
result<-numeric(M)
pi<-3.1415926
result<-sqrt(2)*sin((1:M)*pi*t)
return(result)
}


######
## (ii) cubic b-spline: "poly"
##
B.one<-function(x){
return(x^3/6)
}

##
B.two<-function(x){
result<-(-3*x^3+3*x^2+3*x+1)/6
return(result)
}

##
B.three<-function(x){
result<-(3*x^3-6*x^2+4)/6
return(result)
}

##
B.four<-function(x){
result<-(1-x)^3/6
return(result)
}

##
B.cubic<-function(x){
 result<-0

 if(x>=-3&&x<=-2){
 result<-B.one(x+3)
 }
 
 if(x>=-2&&x<=-1){
 result<-B.two(x+2)
 }

 if(x>=-1&&x<=0){
 result<-B.three(x+1)
 }

 if(x>=0&&x<=1){
 result<-B.four(x)
 }

 return(result)
} 

##
BC.basis<-function(t,knots){
##para:t--evaluation point; knots: equally spaced knots 
M<-length(knots)
delta<-knots[2]-knots[1]
t.v<-(t-knots)/delta

result<-apply(matrix(t.v),1,B.cubic)
return(result)
}

##make it orthonormal
BC.basis.orth<-function(knots,grids){
bs<-apply(matrix(grids),1,BC.basis,knots=knots)
temp<-bs%*%t(bs)*(grids[2]-grids[1])
R<-t(chol(temp))
result<-solve(R, bs)
return(list(result,R))
}

##auxiliary: get R.inv
R.inverse<-function(M,grids=seq(0,1,1/500)){
lmin<-0
lmax<-1 
delta<-(lmax-lmin)/(M-3)
knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
R<-BC.basis.orth(knots,grids)[[2]]
result<-solve(R)
return(result)
}

##evaluate at one points
BC.orth<-function(t,knots,R.inv){
bs.t<-BC.basis(t,knots)
result<-R.inv%*%bs.t
return(result)
}

##orthogonalized cubic B-spline on [0,1] with total M basis (equally spaced knots) 
Basis.Poly<-function(t,M,R.inv){
delta<-1/(M-3)
knots<-c(seq(0,1,delta), 1+delta,1+2*delta)
result<-BC.orth(t,knots,R.inv)
return(result)
}


######
##(iii) natural spline: "ns"
##
NS.orth<-function(grid,df){
 B <- ns(grid, df = df)
 B.orth <- qr.Q(qr(B))  #### orthonormalizing the columns of B 
 return(B.orth)
}

##natural sp;ine evaluate at one points
Basis.NS<-function(t,grid,B.orth){
 ##find the index for t on the grid
   index<-floor(t*length(grid))+1
   result<-B.orth[index,]  
 return(result)
}

######
## (iv) spike basis: three dimensional: "spike"
##
Basis.Spike<-function(t,R.inv){
#a=c(350,220,160)
#b=c(600,190,40)
#c=c(0.7,0.5,0.5)
a<-c(50,70,160)
b<-c(600,300,40)
c<-c(1,1,10)

bs<-matrix(Spike(t,a,b,c))
result<-R.inv%*%bs
return(result)
}



##spike functions: a,b,c: 3 by 1 vectors; a, b: for beta;c: coefficients
Spike<-function(t,a,b,c){
f<-1/beta(a,b)*t^(a-1)*(1-t)^(b-1)
f[3]<-exp(-40*(t-1/2)^2)
result1<-sum(c*f)

pi<-3.1415926
result2<-0.6*sin(4*pi*t)
#result3<-f[3]*sin(32*pi*t)
result3<-f[3]*sin(4*pi*t)
result<-c(result1,result2,result3)
return(result)
}

##get orthogonal transformation matrix for Spike
Spike.orthM<-function(a,b,c,grids){
temp<-apply(matrix(grids),MARGIN=1,Spike,a=a,b=b,c=c)
M<-temp%*%t(temp)*(grids[2]-grids[1])
R<-t(chol(M))
R.inv<-solve(R)
return(R.inv)
}


#######################################II: projection on basis functions
##
Proj.basis.func<-function(psi.v,phi.m,delta){
##para: psi.v: vector of a function evaluated on a fine grid, M by 1
##phi.m: matrix of a basis evaluated on the same fine grid: M by r 
##delta: spacing of the grid
##return: projection (coefficient) of the function on the basis: r by 1
result<-t(psi.v)%*%phi.m*delta
return(result)
}

##
Proj.basis.surf<-function(C.m,phi.m, delta){
##para: C.m: matrix of a surface evaluated on a fine grid^2, M by M
##phi.m: matrix of a basis evaluated on the same fine grid: M by r 
##delta: spacing of the grid
##return: projection (coefficient) of the surface on the basis: r by r matrix
result<-t(phi.m)%*%C.m%*%phi.m*delta^2
return(result)
}

##projection to get inital values for B and lambda
Proj.Ini<-function(covmatrix,M,r,basis.method,grids){
 
 if(basis.method=="poly"){
  R.inv<-R.inverse(M,grids)
  lmin<-0
  lmax<-1 
  delta<-(lmax-lmin)/(M-3)
  knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
  bs<-apply(matrix(grids),1, BC.orth,knots=knots,R.inv=R.inv)
  bs<-t(bs)
 }

 if(basis.method=="ns"){
  bs<-NS.orth(grids,df=M)/sqrt(grids[2]-grids[1])
 }
##
 H<-Proj.basis.surf(covmatrix,bs, delta=grids[2]-grids[1])
 eigenH<-eigen(H,symmetric=TRUE)
 lam.loc<-eigenH$values[1:r]
 B.loc<-eigenH$vectors[,1:r]
 eigen.est2<-t(bs%*%B.loc)
 return(list(B.loc,lam.loc,eigen.est2))
}

###
EigenC<-function(covmatrixd,indext){
##get eigenfunctions and eigenvalues given the covariance surface G(s,t)
##name:EigenC
##para:covmatrixd--estimation of G(s,t);indext--evaluation grid
##result:eigenfunctions and eigenvalues
eigend<-svd(covmatrixd)            ##svd decomposition on the grid of covariance:X=U%*%D%*%t(U)
eigenvd<-eigend$d                  
eigenvd<-eigenvd*(indext[2]-indext[1])  ##eigenvalues
eigenfd<-eigend$u                  
eigenfd<-eigenfd/sqrt(indext[2]-indext[1]) ##eigenfd[,i] is the ith eigenfunction
return(list(eigenfd,eigenvd))
}
###########################################III: Generte eigenfunctions and eigenvalues
###eigen functions: orthonormal; Psi=t(B)%*%Phi, B: M by r; Phi: M by m_i  
Eigenf<-function(t,B,basis.method="sin",R.inv=NULL,grid=NULL){
##para: t--evaluation point; B: coefficient: M by r; 
##where M--dimension of basis to use; r: dimension of process; M >=r
##basis.method: which basis to use; one of "sin", "poly", "spike", "ns"
##return: eigenfunctions evaluated  at t; r by 1 
 M<-nrow(B)
 r<-ncol(B)

result<-numeric(r)
 if(basis.method=="sin"){
 phi.c<-Basis.Sin(t,M)
 }

 if(basis.method=="poly"){
 phi.c<-Basis.Poly(t,M,R.inv)
 }
 
 if(basis.method=="spike"){
 phi.c<-Basis.Spike(t,R.inv)[1:M]
 }
 
 if(basis.method=="ns"){
 phi.c<-Basis.NS(t,grid,R.inv)
 }
 result<-t(B)%*%matrix(phi.c)
 return(result)
}


##eigenvalues (>=0)
Eigenv<-function(r, method="algebra",alpha, beta=1,r1=5,r2=10, gamma=1){
##para: r--number of nonzero eigenvalues
##method: method to generate: one of "algebra", "geometric", "hyperbolic", "hybrid"
##alpha, beta
##return: first r eigenvalues; r by 1

 result<-numeric(r)
 if(method=="algebra")
   temp<-(1:r)^(-alpha)

 if(method=="geometric")
   temp<-(alpha)^(1:r)

  if(method=="hyperbolic")
   temp<-(1-alpha/(1:r))^(-beta)
 
  if(method=="hybrid"){
   temp1<-(1:r1)^(-alpha)
   temp2<-temp1[r1]*exp(-gamma*(1:r2))    
   temp<-c(temp1,temp2)
  }

##rescale such that the first (also the largest) eigenvalue is  one
  result<-temp/temp[1]
  return(result)
}


######################################## IV: Generate Curves
###(i)random curve with mu(t)=0; under normality assumption, i.e., ksi i.i.d N(0,1)
RanCur<-function(t.v,r,eigenf.method="sin",B,eigenv.method="algebra",alpha, beta=1,r1=5,r2=10,gamma=1,R.inv=NULL,grid=NULL){
##t.v--vector of evaluation points 
##r--number of nonzero eigenvalues
##eigenf.method: the way to generate eigenfunctions; one of "sin", "poly",
##eigenv.method: the way to generate eigenvalues; one of "algebra", "geometric" and "hyperbolic"
##return: random curve evaluated at t.v: length(t.v) by 1; and pca scores: r by 1

result<-numeric(length(t.v))

eigen.f<-apply(matrix(t.v),MARGIN=1,Eigenf,basis.method=eigenf.method,B=B,R.inv=R.inv,grid=grid)

eigen.v<-matrix(Eigenv(r,eigenv.method,alpha,beta,r1,r2,gamma),r,length(t.v))

ksi.v<-rnorm(r)
ksi.m<-matrix(ksi.v,r,length(t.v))

temp<-sqrt(eigen.v)*eigen.f*ksi.m
result<-apply(temp,2,sum)

return(list(result,ksi.v))
}


############Part II: simulation
##measurement points: N from uniform [nmin,nmax];L from a beta(a,b) distribution;
MeasL<-function(nmin,nmax,a=1,b=1){
##para:N from uniform [nmin,nmax];L from a beta(a,b) distribution;
##return: measurement points;
N<-sample(nmin:nmax,1)
result<-sort(rbeta(N,a,b))
return(result)
}


##generate one sparsely observed realization with measurement error; 
##also the corresponding random curve evaluated on a fine grid (omit this)
GenObs<-function(r,eigenf.method="sin",B,eigenv.method="algebra",alpha,beta=1,nmin,nmax,a=1,b=1,grid,sig,R.inv=NULL,grid.ns=NULL,r1=5,r2=10,gamma=1){
 L.v<-MeasL(nmin,nmax,a,b)
# t.v<-c(L.v,grid)
 t.v<-L.v
 e.v<-rnorm(length(L.v))*sig

 result.all<-RanCur(t.v,r,eigenf.method,B, eigenv.method,alpha, beta,r1,r2,gamma,R.inv,grid.ns)

 Obs.all<-result.all[[1]]
 Ksi.all<-result.all[[2]]

 Obs.v<-Obs.all[1:length(L.v)]+e.v  
 result1<-cbind(Obs.v,L.v)

# Obs.grid<-Obs.all[-(1:length(L.v))] 
# result2<-cbind(Obs.grid,grid) 
# return(list(result1,result2,Ksi.all))  
  return(list(result1,Ksi.all))  

}


###
GenObs.T<-function(r,eigenf.method="sin",B,eigenv.method="algebra",alpha,beta=1,nmin,nmax,a=1,b=1,grid,sig,R.inv=NULL,df=4,factor=1.413,r1=5,r2=10,gamma=1){
 L.v<-MeasL(nmin,nmax,a,b)
#t.v<-c(L.v,grid)
 t.v<-L.v
 e.v<-rt(length(L.v),df)*sig/factor  ### noise generation from t distribution

 result.all<-RanCur(t.v, r, eigenf.method, B, eigenv.method,alpha, beta,r1,r2,gamma,R.inv)

 Obs.all<-result.all[[1]]
 Ksi.all<-result.all[[2]]

 Obs.v<-Obs.all[1:length(L.v)] + e.v  ### noise addition to the true sample curves
 result1<-cbind(Obs.v,L.v)

# Obs.grid<-Obs.all[-(1:length(L.v))] 
# result2<-cbind(Obs.grid,grid) 
# return(list(result1,result2,Ksi.all))  
  return(list(result1,Ksi.all)) 

}
#################################
### For t-distribution
## d.f. 3 : s.d. = 1.725
## d.f. 4 : s.d. = 1.413
## d.f. 5 : s.d. = 1.291
## d.f. 6 : s.d. = 1.224
## d.f. 7 : s.d. = 1.183
## d.f. 8 : s.d. = 1.155
## d.f. 9 : s.d. = 1.134
## d.f. 10: s.d. = 1.118
## d.f. 12: s.d. = 1.096
## d.f. 15: s.d. = 1.074










