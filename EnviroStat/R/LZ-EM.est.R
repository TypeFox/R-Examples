

staircase.EM =function(data,p=1,block=NULL,covariate=NULL,B0=NULL,init=NULL,a=2,r=.5,verbose=FALSE,
     maxit=20,tol=1e-6) {

#########################
# Model:
#      data  ~   MVN ( z x beta , kronecker(I, Sigma) )
#      beta  ~   MVN (Beta0 , kronecker(Finv , Sigma )
#      Sigma ~   GIW (Theta , delta's )
#  
#      Theta is a collection of hyperparamaters including Xi0, Omega, Lambda, Hinv
# 
# Input:
#  data: data matrix, grouped by blocks each with stations having the same number of 
#     missing observations. The blocks are organized in order of decreasing number 
#     of missing observations, ie. block 1 has more missing observations than block 2.
#  Default structure: 
#     - Each column represent data from a station; rows are for time
#     - Blocks are decided based on the number of missing observations

# Optional input:
#  p: number of pollutants measured at each stations. 
#  (first p columns of y are for p pollutants from station 1, block 1).
#   
#  block: a vector indicating the number of stations in each block - from 1 to K
#  covariate: design matrix for covariates created with "model.matrix" with "as.factor" 
#  B0: Provided if the hyperparameter beta_0 (B0) is known and not estimated
#  init: Initial values for the hyperparameters; output of this function can be used for that
#  a,r : When p=1, the type-II MLE's for delta's are not available. Delta's are assumed to follow
#        a gamma distribution with parameters (a,r)
#  verbose: flag for writing out the results at each iteration
#  maxit: the default maximum number of iteration
#  tol: The convergence level.
# 
# Output:
#  Delta:  The estimated degrees freedom for each of the blocks (list)
#  Omega:  The estimated covariance matrix between pollutants
#  Lambda: The estimated conditional covariance matrix between stations in each block given
#           data at stations in higher blocks (less missing data) - (list)
#  Xi0:    The estimated slopes of regression between stations in each blocks and those in higher
#            blocks (list). Note that tau_0i = kronecker(Xi0, diag(p)) - same across stations
#            for each pollutants.
#  Beta0:  Coefficients - assumed to be the same across stations for each pollutant
#  Finv:   Scale associated with beta_0
#  Hinv :  The estimated hyperparameters (list) - inverse of H_j
#  Psi:    The estimated (marginal) covariance matrix between stations     
#  block:  From input
#  data:   From input
#  covariate: From input
#
#  Lambda.1K: The inverse Bartlett decomposition (Le+Zidek 2006, Section 15.2, pages 302-303)
#   
##########################
# Details on estimation given in Le and Zidek (2006 - Statistical Analysis of Environmental Space-Time Process, Springer)
##########################
 

# Some useful functions for the estimation procedure

# A useful identity for obtaining inverse of (A + B'CB) using Ainv and Cinv via
# (A + BCB')^{-1}= Ainv - Ainv*B (Cinv+B'*Ainv*B)^{-1} B'*Ainv
invert=function(Ainv,B,Cinv) {
  tmp=solve(Cinv+crossprod(B,crossprod(Ainv,B)),crossprod(B,Ainv))
  tmp=crossprod(t(crossprod(Ainv,B)),tmp)
  return(Ainv-tmp)
}
# Above identity with A = Identity
invAeI=function(B,Cinv) {
  tmp=crossprod(t(B),solve(Cinv+crossprod(B),t(B)))
  return(diag(dim(tmp)[1])-tmp)
}

# To obtain determinant of a large non-negative matrix
# through the Bartletts decomposition

det.BD=function(x,max.size=20,logarithm=TRUE) {
  val=0
  dm=dim(x)[1]
  n=ceiling(dm/max.size)
  i2=floor(seq(0,dm,length=n+1))
  i1=(i2+1)[1:n]
  i2=i2[1:n+1]
  for (i in n:1) {
    x22=x[i1[i]:i2[i],i1[i]:i2[i],drop=FALSE]
    if (i>1) {
      x21=x[i1[i]:i2[i],1:i2[i-1],drop=FALSE]
      x11=x[1:i2[i-1],1:i2[i-1],drop=FALSE]
      x=x11-crossprod(x21,solve(x22,x21))
    }
    val=val+determinant(x22)[[1]][1]
  }
if (logarithm) return(val)
else return(exp(val))
}

# evaluating log density of matric-t distribution
lden.matrict = function(x,mu=0,A,B,delta) 
{
#A :n x n; B :m x m ; x :n x m

d =delta
n =dim(x)[1]
m =dim(x)[2]
p = n+m
lK =  -n*m*log(d*pi^2)/2 + sum(lgamma((d+p - (1:p))/2)) + p*(p-1)*log(pi)/4 -
     sum(lgamma((d+n - (1:n))/2))-n*(n-1)*log(pi)/4 -
     sum(lgamma((d+m - (1:m))/2))-m*(m-1)*log(pi)/4
#
t1 = determinant(A)[[1]][1]
t2 = determinant(B)[[1]][1]
t3 = determinant(diag(n)+solve(A)%*%(x-mu)%*%t((x-mu)%*%solve(B))/d)[[1]][1]
lK = lK - m*t1/2 - n*t2/2 - (d+p-1)*t3/2
return(lK)
}

########################################################################
#                 Program begins here
########################################################################
# Create lists needed
Ft=list()
Fhat=list()
W=list()
Bt=list()
EGA=list()
EGAI = list()
ESI=list()
LAM=list()
LAM.temp=list()
PSIt=list()
PSI=list()
PHI=list()
SSQ=list()
Hinv=list()
lambdastar=list()
Ht=list()
T0=list()
T0t=list()
Bhat=list()
EBS=list()
XI0=list()
XI0.temp=list()
EGTH=list()
OMEGA=NULL
OMEGA.temp=NULL

y = data
covariate = as.matrix(covariate)  #**** for only one covariate provided as vector, eg. overall mean
z = covariate
# Processing data provided and 
#    compute dimensions

if (class(y)=="data.frame") y=as.matrix(y)
ii=apply(!is.na(y),1,sum)!=0
y=y[ii,]
n=dim(y)[1]
if (!is.null(z) & length(z)>1) {
if (length(z)==length(ii))
  z=matrix(z[ii],ncol=1)
else if (dim(z)[1]==length(ii))
  z=z[ii,]
}
g1=table(apply(is.na(y),2,sum))/p
g1 = g1[length(g1):1]
if (is.null(block)) block=g1
g = block
if (sum(g1)!=sum(g))
stop(paste("You specified",sum(g),"sites but data appears to have",sum(g1),"sites"))
if (!all(cumsum(g1)%in%cumsum(g)))
stop("All sites within a cluster must have the same pattern of missing values")
m0=as.numeric(names(g1))
m=g
ii=match(cumsum(g1),cumsum(g))
ii=cbind(c(1,ii[-length(ii)]+1),ii)
for (i in 1:length(m0))
m[ii[i,1]:ii[i,2]]=m0[i]
g=c(g)
d=g*p + 10
if (p==1) d=g*p + a/r
K=length(d)
dm=cbind(c(0,p*cumsum(g[-K]))+1 ,p*cumsum(g))
row.names(dm) = names(g1)
names(block) = c(1:K)
ELG=numeric(K)



# Setting Z and evaluating Bhat and Fhat=t(Zj)%*%Zj 
#   as well as initializing B0 and FF, delta

if (is.null(z)) {
Finv=NULL
B0 = NULL
epst = y
}
else {
b0 = is.null(B0)
Z = z
# Initialize B0
if (b0) B0 = NULL
FF= matrix(0,dim(Z)[2],dim(Z)[2])
for (i in K:1) {
FF=FF+g[i]*(Fhat[[i]]=crossprod(Z[(m[i]+1):n,]))
FhatInv = ginv(Fhat[[i]])
# Get Bhat[i,...,k]; see Le+Zidek (2006) Section 10.4
Bhat[[i]] = crossprod(FhatInv,crossprod(Z[(m[i]+1):n,],y[(m[i]+1):n,dm[i,1]:dm[K,2]]))
# Bhat_0^j
if (b0) { 
 if (dim(Z)[2]==1) B0 = c(Bhat[[i]][,1:(g[i]*p)],B0)
   else  B0 = cbind(Bhat[[i]][,1:(g[i]*p)],B0)   
      }
}

q = dim(Bhat[[K]])[1]

if (b0) {     # *** to cover the case when B0 is provided and doesn't need to be estimated
    if (dim(Z)[2]> 1) B0 = matrix(apply(B0,1,mean),nrow=q,ncol=dm[K,2])
    else B0 = matrix(B0,nrow=q,ncol=dm[K,2])
}

if (dim(Z)[2]!=dim(B0)[1])
 stop("B0 has dimensions inconsistent with Z")
priortrend=crossprod(t(Z),B0)
epst = y-priortrend
FF=FF/sum(g)
Finv=solve(FF)
if (b0) R = matrix(rep(diag(p),sum(g)),nrow=p)
}

# Initialize LAMDA - dimension g_i x g_i
for (i in K:1) {

x2=NULL
i1=matrix(1:(g[i]*p),nrow=p)
x1=epst[(m[i]+1):n,dm[i,1]:dm[i,2],drop=FALSE]
for(j in 1:dim(i1)[2])
x2=cbind(x2,c(x1[,i1[,j],drop=FALSE]))
LAM[[i]]=var(x2,na.rm=T)
LAM[[i]]=d[i]*LAM[[i]]
}

# Initialize OMEGA : dimension p x p 
i1=matrix(1:(sum(g)*p),nrow=p)
for (i in 1:dim(i1)[1])
OMEGA=cbind(OMEGA,c(epst[,i1[i,]]))
OMEGA=var(OMEGA,na.rm=T)
OMEGA = OMEGA/max(diag(OMEGA))

# Initialize XI : dimension (g_{i+1} + ... + g_k) x g_i

dm1 = dm/p
dm1[,1]=(dm[,1]-1)/p+1
if (K>1) for (i in (K-1):1) {
x2=NULL
x1=epst[(m[i]+1):n,dm[i,1]:dm[K,2],drop=FALSE]
i1=matrix(1:dim(x1)[2],nrow=p)

# put all columns in each stations into one and blocks i to k together

for (j in 1:dim(i1)[2]) x2=cbind(x2,c(x1[,i1[,j],drop=FALSE])) 
i2 = dm1[i,2] - dm1[i,1] + 1   # number of stations in block i 
i3 = dm1[K,2] - dm1[i,1] + 1   # number of stations in blocks i to k 

XI0[[i]]=solve(crossprod(x2[,(i2+1):i3],x2[,(i2+1):i3]),
        crossprod(x2[,(i2+1):i3],x2[,1:i2]))

# Obtain tau_o

T0[[i]]=kronecker(XI0[[i]],diag(p)) 
}

# check if initial values are supplied

if (!is.null(init)) {
 nm1=names(init)
 if ("Delta" %in% nm1) d=init$Delta
 if ("Omega" %in% nm1) OMEGA=init$Omega
 if ("Lambda" %in% nm1) LAM=init$Lambda
 if ("Xi0" %in% nm1) XI0=init$Xi0
 if (!is.null(z)) {
   if ("Beta0" %in% nm1){ B0=matrix(init$Beta0,ncol=p*sum(g))
                         priortrend=crossprod(t(Z),B0)
                         epst=y-priortrend
                       }
   if ("Finv" %in% nm1) { Finv=init$Finv; FF=solve(Finv) }
  }
}

# Compute marginal density with initial values - Le+Zidek(2006) Section 10.3, pages 159-161
MD1=0
dd0 = d[K]-g[K]*p+1
if (is.null(z)) PHI[[K]]=diag(n)
else PHI[[K]]=invAeI(Z,FF)
A = solve(PHI[[K]])
B=kronecker(LAM[[K]]/(dd0),OMEGA)
x1=epst[(m[K]+1):n,dm[K,1]:dm[K,2],drop=FALSE]
MD1 = MD1 + lden.matrict(x1,,A,B,dd0) 

lambdastar[[K]] = LAM[[K]]  
if (K >1) for (j in (K-1):1) {
 x1=epst[(m[j]+1):n,dm[j,1]:dm[j,2],drop=FALSE] 
 x2 = epst[(m[j]+1):n,dm[j+1,1]:dm[K,2],drop=FALSE] 
 dd0 = d[j]-g[j]*p+1 
 P1=crossprod(lambdastar[[j+1]],XI0[[j]])
 P2=crossprod(XI0[[j]],P1)
 lambdastar[[j]] = rbind(cbind(LAM[[j]]+P2,t(P1)),cbind(P1,lambdastar[[j+1]]))
 Hinv[[j]]=kronecker(lambdastar[[j+1]],OMEGA)
 if(is.null(z)) PHI[[j]]= invAeI(x2,Hinv[[j]])
  else { 
 Ainv=invAeI(Z[(m[j]+1):n,],FF)
 PHI[[j]]=invert(Ainv,x2,Hinv[[j]])
  }
A = solve(PHI[[j]])
B = kronecker(LAM[[j]]/(dd0),OMEGA)
mu = x2 %*% T0[[j]]
# if (verbose) cat(paste("mu ",mu[1,1:3],"\n")) 
MD1 = MD1 + lden.matrict(x1,mu,A,B,dd0)
}
if (p==1) { dden = sum(log(dgamma(d,shape=a,scale=1/r)))
            if (dden < -1000) dden = log(1/10e10) 
            MD1= MD1+ dden
           }
if (verbose) cat(paste("Log likelihood = ",round(MD1,5)," at iteration 0 \n"))

#######################################################
#
# Iterating between the EM steps
#
######################################################

llike = MD1 
for (iter in 1:maxit) {
if (verbose) cat(paste("\nIteration",iter,"\n"))

#####################################################
# E-STEP: Evaluating the expectations given the current estimates 
# Le+Zidek (2006) Section 10.5, pages 165-167
#####################################################

# Starting with the last block (Block K: Full data stations - no T0 or H required)
# Extract data
x1=epst[,dm[K,1]:dm[K,2],drop=FALSE]   

# Evaluate the expectations given para. estimates at previous iter
lambdastar[[K]] = LAM[[K]]  # this is used to evaluate Hinv for other blocks
PSI[[K]]=LAM[[K]]/(d[K]-g[K]*p-1)
PHI[[K]]=diag(n)  
if (is.null(Finv)) SSQ[[K]]=crossprod(x1)
else {
  PHI[[K]]=invAeI(Z,FF) #A22
  SSQ[[K]]=crossprod(x1,crossprod(PHI[[K]],x1))
}
PSIt[[K]]=kronecker(LAM[[K]],OMEGA)+SSQ[[K]]
ESI[[K]]=EGAI[[K]]=(d[K]+n-m[K])*solve(PSIt[[K]])

ELG[K]=det.BD(PSIt[[K]])-g[K]*p*log(2)-sum(digamma((d[K]+n-m[K]+1-1:(g[K]*p))/2))

# Evaluate components of tilde(beta_K) and E(beta Sig^{-1}), 
#  and E(beta Sig^{-1} beta) - corresponding to block K.
# Note: These are only needed if covariate(s) exists
if (!is.null(Finv)) {
Ft[[K]]=Fhat[[K]]+FF  # Eqn(12)
FtInv = ginv(Ft[[K]])
W[[K]]=crossprod(FtInv, Fhat[[K]]) 
Bt[[K]]=crossprod(t(W[[K]]),Bhat[[K]])+
       crossprod(t(diag(q)-W[[K]]),B0[,dm[K,1]:dm[K,2],drop=FALSE])
EBS[[K]]=crossprod(t(Bt[[K]]),EGAI[[K]]) # (11) the first component
# computing (17) the first component  
EBSB=crossprod(t(EBS[[K]]),t(Bt[[K]]))+g[K]*p*solve(Ft[[K]])   
}


# If number of blocks K > 1,
# Continue with other blocks -Block 1 to (K-1): The rest of the staircase
#
if (K>1) for ( i in (K-1):1) {

# Evaluate Hinv_i (assuming the IW pattern)
P1=crossprod(lambdastar[[i+1]],XI0[[i]]) 
P2=crossprod(XI0[[i]],P1) 
lambdastar[[i]] = rbind(cbind(LAM[[i]]+P2,t(P1)),cbind(P1,lambdastar[[i+1]]))
Hinv[[i]]=kronecker(lambdastar[[i+1]],OMEGA)


# Evaluate Psi
P1=crossprod(PSI[[i+1]],XI0[[i]])
P2=crossprod(XI0[[i]],P1)
trpsih = sum(diag(crossprod(PSI[[i+1]],solve(lambdastar[[i+1]]))))*p
P3=LAM[[i]]*(1+trpsih)/(d[i]-g[i]*p-1)
PSI[[i]]=rbind(cbind(P3+P2,t(P1)),cbind(P1,PSI[[i+1]]))

# Evaluate tilde(tau), tilde(H_i), and component A22
# Extract data from block i (x1) and  blocks i+1 to K (x2)
x1=epst[(m[i]+1):n,dm[i,1]:dm[i,2],drop=FALSE]  
x2=epst[(m[i]+1):n,dm[i+1,1]:dm[K,2],drop=FALSE]

# if no covariates 
if (is.null(Finv)) {
Ht[[i]]=solve(Hinv[[i]]+crossprod(x2))
T0t[[i]]=crossprod(Ht[[i]],(crossprod(Hinv[[i]],T0[[i]])+
        crossprod(x2,x1)))
PHI[[i]]=invAeI(x2,Hinv[[i]])
}
# if there are covariates
else {
Ainv=invAeI(Z[(m[i]+1):n,],FF)
Ht[[i]]=solve(Hinv[[i]]+crossprod(x2,crossprod(Ainv,x2)))
T0t[[i]]=crossprod(Ht[[i]],(crossprod(Hinv[[i]],T0[[i]])+
        crossprod(x2,crossprod(Ainv,x1)))) # (6)
PHI[[i]]=invert(Ainv,x2,Hinv[[i]])
}
x3=x1-crossprod(t(x2),T0[[i]])
SSQ[[i]]=crossprod(x3,crossprod(PHI[[i]],x3))
PSIt[[i]]=kronecker(LAM[[i]],OMEGA)+SSQ[[i]]
EGAI[[i]]=(d[i]+n-m[i])*solve(PSIt[[i]])

# Evaluate E[SIGMA^{-1}{i,i}] in Eqn 23
tmp=crossprod(t(T0t[[i]]),EGAI[[i]]) # 1,2 and 2,1 pieces of E[SIGMA^{-1}{i,i}]
tmp1=crossprod(t(tmp),t(T0t[[i]]))+g[i]*p*Ht[[i]]+ESI[[i+1]] # 2,2 piece of  E[SIGMA{i,i}]
ESI[[i]]=rbind(cbind(EGAI[[i]],-t(tmp)),cbind(-tmp,tmp1))

# E(Gamma^{-1} tau')Hinv 
EGTH[[i]] = crossprod(tmp,Hinv[[i]])


# Evaluate E(log(|Gamma_i|)
ELG[i]=det.BD(PSIt[[i]])-g[i]*p*log(2)-sum(digamma((d[i]+n-m[i]+1-1:(g[i]*p))/2))

# Adding components to E(beta Sig^{-1} beta)
if (!is.null(Finv)) {
Ft[[i]]=Fhat[[i]]+FF # (12)
FtInv = ginv(Ft[[i]])
W[[i]]=crossprod(FtInv, Fhat[[i]]) # between (12) & (13) 
Bt[[i]]=crossprod(t(W[[i]]),Bhat[[i]])+
      crossprod(t(diag(q)-W[[i]]),B0[,dm[i,1]:dm[K,2],drop=FALSE])

tmp=Bt[[i]][,1:(g[i]*p),drop=FALSE]-
  crossprod(t(Bt[[i]][, (g[i]*p+1):(dim(Bt[[i]])[2]),drop=FALSE]),T0t[[i]])
tmp0=crossprod(t(tmp),EGAI[[i]])
tmp1=crossprod(t(tmp0),t(T0t[[i]]))
tmp2=g[i]*p*crossprod(t(Bt[[i]][,(g[i]*p+1):(dim(Bt[[i]])[2]),drop=FALSE]),Ht[[i]])
# Augmenting E(BSigma^{-1}) - Eqn 11
EBS[[i]]=cbind(tmp0, EBS[[i+1]]-tmp1+tmp2)  
# Adding relevant components to E(BSigma^{-1} B)
EBSB=EBSB+g[i]*p*solve(Ft[[i]])+
  crossprod(t(tmp2),t(Bt[[i]][,(g[i]*p+1):(dim(Bt[[i]])[2]),drop=FALSE]))
EBSB=EBSB+crossprod(t(tmp0),t(tmp))
}
}

#############################################################
#
# M-STEP: Maximize the objective function to get new estimates of parameters
# See Section 10.6 - Maximization Equations, Le+Zidek(2006)
#############################################################
#  Need iteration between Omega + Lambda and delta
# 
LAM.temp = LAM
OMEGA.temp = OMEGA
d.temp = d
for (iter.sub1 in 1:10) {

# New estimates of Lambda[[i]] for i = 1 to K
for (i in K:1) {
C=matrix(0,g[i],g[i])
ii=matrix(1:(p*g[i]),ncol=g[i])
for (i0 in 1:g[i]) for (j0 in 1:g[i])
	C[i0,j0]=sum(OMEGA.temp*(EGAI[[i]][ii[,i0],ii[,j0],drop=FALSE]))
A = matrix(d.temp[i]*p*solve(C),g[i],g[i]) 
LAM.temp[[i]]= (A +t(A))/2
}

# New estimate of Omega

C=array(0, c(p, p, K))
for(i in 1:K) {
ii=matrix(1:(p*g[i]),byrow=T,ncol=p)
for(i0 in 1:p) for(j0 in 1:p) 
C[i0,j0,i]=sum(LAM.temp[[i]]*(EGAI[[i]][ii[,i0],ii[,j0],drop=FALSE]))
}
A =sum(g*d.temp)*solve(apply(C,c(1,2),sum))
OMEGA.temp = (A +t(A))/2

# New estimates of Delta's

if (p > 1) {
a = 1
r = 0
}
for (i in 1:K) {
for (l in 1:30) {
f0=p*det.BD(LAM.temp[[i]])+g[i]*det.BD(OMEGA.temp)-ELG[i]-p*g[i]*log(2)-
	sum(digamma((d.temp[i]+1-1:(g[i]*p))/2))+(a-1)/d.temp[i]-r
f2=sum(-trigamma((d.temp[i]+1-1:(g[i]*p))/2))/2-(a-1)/d.temp[i]^2
if (abs(f0)<tol/100) break
d.temp[i]=max(d.temp[i]-f0/f2,p*g[i]+2)

}
 }
}


# New estimate of XI0  

  if (K>1) for (i in (K-1):1) 
  {
  ntmp = dm1[K,2] - dm1[i,2]
  G=matrix(0, g[i]*ntmp, g[i]*ntmp)
  i0=matrix(1:(g[i]*p),byrow=T,ncol=p)
  j0=matrix(1:(ntmp*p),byrow=T,ncol=p)
  for (ii in 1:p) for (jj in 1:p) {
 	B=EGAI[[i]][i0[,ii],i0[,jj],drop=FALSE]
  C=Hinv[[i]][j0[,ii],j0[,jj],drop=FALSE]
 	G=G+kronecker(B,C)
 	}
 D=rep(0,ntmp*g[i])
 for(ii in 1:p) 
 	D=D+c(t(EGTH[[i]][i0[,ii],j0[,ii],drop=FALSE]))
 M=solve(G,D)
 XI0.temp[[i]]=matrix(M,ncol=g[i])
 }


# Update B0 and FF
# 

if (!is.null(Finv)) {
if (b0) {
# get exchangeable structure between stations - ie. different for
#   each of the components in the multivariate response
beta.tmp = solve(crossprod(t(R),crossprod(ESI[[1]],t(R))), crossprod(t(R),t(EBS[[1]])))
B0.tmp = crossprod(beta.tmp,R)
}
else B0.tmp = B0  # ***  to cover the case when B0 is provided
tmp0=crossprod(t(EBS[[1]]),t(B0.tmp))
tmp1=crossprod(t(B0.tmp),crossprod(ESI[[1]],t(B0.tmp)))
Finv.tmp=(EBSB-tmp0-t(tmp0)+tmp1)/dm[K,2]
FF.tmp=solve(Finv.tmp)
}

# Replacing the current estimates with new ones for next iteration

OMEGA = OMEGA.temp
d = d.temp
LAM = LAM.temp
if (K > 1) {
XI0 = XI0.temp
# Obtain tau_o
for (i in (K-1):1)
T0[[i]]=kronecker(XI0[[i]],diag(p)) 
}

if (!is.null(z)) {
if (b0) B0 = B0.tmp
priortrend=crossprod(t(Z),B0)
epst=y-priortrend
Finv = Finv.tmp
FF = FF.tmp
}

# Compute marginal density after each iteration
MD1=0
dd0 = d[K]-g[K]*p+1
if (is.null(z)) PHI[[K]]=diag(n)
else PHI[[K]]=invAeI(Z,FF)
A = solve(PHI[[K]])
B=kronecker(LAM[[K]]/(dd0),OMEGA)
x1=epst[(m[K]+1):n,dm[K,1]:dm[K,2],drop=FALSE]
MD1 = MD1 + lden.matrict(x1,,A,B,dd0)

lambdastar[[K]] = LAM[[K]] 
if (K >1) for (j in (K-1):1) {
 x1=epst[(m[j]+1):n,dm[j,1]:dm[j,2],drop=FALSE]
 x2 = epst[(m[j]+1):n,dm[j+1,1]:dm[K,2],drop=FALSE]
 dd0 = d[j]-g[j]*p+1
 P1=crossprod(lambdastar[[j+1]],XI0[[j]])
 P2=crossprod(XI0[[j]],P1)
 lambdastar[[j]] = rbind(cbind(LAM[[j]]+P2,t(P1)),cbind(P1,lambdastar[[j+1]]))
 Hinv[[j]]=kronecker(lambdastar[[j+1]],OMEGA)
 if(is.null(z)) PHI[[j]]= invAeI(x2,Hinv[[j]])
  else { 
    Ainv=invAeI(Z[(m[j]+1):n,],FF)
    PHI[[j]]=invert(Ainv,x2,Hinv[[j]])
  }
A = solve(PHI[[j]])
B = kronecker(LAM[[j]]/(dd0),OMEGA)
mu = x2 %*% T0[[j]]
MD1 = MD1 + lden.matrict(x1,mu,A,B,dd0)
}
if (p==1) { dden = sum(log(dgamma(d,shape=a,scale=1/r)))
            if (dden < -1000) dden = log(1/10e10) 
            MD1= MD1+ dden
           }

# Print interm results and check for convergence

if (verbose) for (i in 1:K) cat(paste(" Delta",i,"=",d[i],"\n"))
if (verbose) cat(paste("Log likelihood = ",round(MD1,5),"Relative change = ",abs((MD1-llike)/llike),"\n"))
if (abs((MD1-llike)/llike)<tol) break
llike = MD1

if (iter == maxit) {
	cat(paste("Maximum (",maxit,") iterations exceeded without convergence\n"))
}
}
if (verbose) cat("\n\n")
names(d) = NULL
del = list()
for (i in 1:K) del[[i]] = d[i]
B0 = as.matrix(B0) #*** to ensure B0 is passed as matrix
obj=list(Delta=del,Omega=OMEGA,Lambda=LAM,Xi0=XI0,Beta0=B0,Finv=Finv,Psi=PSI,
                     Hinv=Hinv,data=data, block= block, covariate=covariate)
return(obj)
}



staircase.hyper.est <- function(emfit,covfit,u,p,g,d0=NULL)
{
# This function combines the results from the "staircase.EM" fit and the
# SG method to estimate the hyperparameters associated with the ungauged sites

#Input
# emfit : Output from the staircase.EM fit
# covfit : The covariance matrix between all locations (with new locations 
#          at the beginning). This is an output from the SG fitting
# u: number of new locations
# p: dimension of the multivariate response
# g: number of stations
# d0 (Optional): The degrees of freedom for the new locations (ungauged block)
#Output
#      Delta.0: The degree of freedoms for the new locations 
#               = d0 if given (must be > u*p+2)
#               else
#			= mean(emfit$delta) if > u*p+2
#          		= u*p+ min(emfit$delta) otherwise
#      Lambda.0 : Conditional variance between new locations given the gauged stations
#      Xi0.0 : the regression slope (Note: Tau_0i = Kronecker(Xi0 , diag(p))
#      H.0 : The variance matrix for the rows of Tau^[u]
#  Also all components of the output of the staircase.EM fit (for blocks 1-K).

if (!is.null(d0)) {
 if (d0 < u*p +2) stop(paste("The provided degrees of freedom (",d0,") is smaller than", u*p+2))
 }
omega = emfit$Omega
K= length(emfit$block)
g = sum(emfit$block)
delta = NULL
for (i in 1:K) 
 delta = c(delta, emfit$Delta[[i]])

lambdastar = list()
lambdastar[[K]] = emfit$Lambda[[K]]  
if (K >1) for (j in (K-1):1) {
 P1=crossprod(lambdastar[[j+1]],emfit$Xi0[[j]])
 P2=crossprod(emfit$Xi0[[j]],P1)
 lambdastar[[j]] = rbind(cbind(emfit$Lambda[[j]]+P2,t(P1)),cbind(P1,lambdastar[[j+1]]))
 }

psi = emfit$Psi[[1]]
if (is.null(d0)) {
delta.u = u*p + min(delta)
if (mean(delta) > delta.u) delta.u = mean(delta)
 }
else delta.u = d0
h.u = kronecker(solve(lambdastar[[1]]),solve(omega))
xi.u  = solve(covfit[(u+1):(u+g),(u+1):(u+g)]) %*% covfit[(u+1):(u+g),1:u]

tr = sum(diag(kronecker(psi,omega)%*%h.u))
lambda.u = (delta.u - u*p -1) * (covfit[1:u,1:u] - t(xi.u)%*%
        covfit[(u+1):(u+g),(u+1):(u+g)] %*% xi.u )/ (1+ tr)


obj = list(Delta.0 = delta.u, Delta = emfit$Delta, Lambda.0 = lambda.u, Lambda=emfit$Lambda,
  Omega= omega, Xi0=emfit$Xi0, Xi0.0 = xi.u, Beta0=emfit$Beta0,Finv=emfit$Finv, Psi=emfit$Psi,
  H.0= h.u,Hinv=emfit$Hinv, data=emfit$data, block= emfit$block, covariate=emfit$covariate)

return(obj)
}




