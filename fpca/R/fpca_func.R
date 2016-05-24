##################################
######02-06-08: R package "fpca"
######R functions for MLE of FPCA to longitudinal data
##name: fpca.mle
##prupose: use Newton's method to fit MLE for FPCA, as well as do model selection by 
##         the approximate CV score
fpca.mle<-function(data.m, M.set,r.set,ini.method="EM", basis.method="bs",sl.v=rep(0.5,10),max.step=50,
grid.l=seq(0,1,0.01),grids=seq(0,1,0.002)){
##para:data.m: the data matrix with three columns: column 1: subject ID, column 2: observation, column 3: measurement time
##M.set--# of basis functions; r.set: dimension of the process
##ini.method: initial method for Newton, one of "EM","loc"; 
##basis.method: basis functions to use for Newton, one of "bs" (cubic Bsplines), "ns" (nature splines)
##sl.v: shrinkage steps for Newton, e.g., sl.v<-c(0.5,0.5,0.5) means the first three steps are 0.5
##max.step: max number of iterations of Newton
##grid.l: grids for loc; grids: denser grids for EM/Newton;
##return: a list of (i) the selected model; (ii) the corresponding eigenfunctions; (iii)eigenvalues
##        (iv) error variance; (v) fitted mean (by local linear); (vi) the grid where the evaluation takes place

##(i)set some parameters in fpca.fit
tol<-1e-3                                        #tolerance level to determine convergence in Newton 
cond.tol<-1e+10                                  #tolerance to determine singularity in Newton 
if(length(sl.v)<max.step){
sl.v<-c(sl.v,rep(1,max.step-length(sl.v)))
}else{
sl.v<-sl.v[1:max.step]
}

if(basis.method=="bs"){
basis.method<-"poly"
}

M.self<-20
basis.EM<-"poly"                                  ##basis for EM

if(ini.method=="EM"){
no.EM<-FALSE
}else{
no.EM<-TRUE
}

iter.EM.num<-50                                   ##maximum number of iterations of EM 
sig.EM<-1                                         ##initial value of sig in EM 


nmax<-max(table(data.m[,1]))   ##maximum number of measurements per subject  
L1<-min(as.numeric(data.m[,3]))              ##range of the time 
L2<-max(as.numeric(data.m[,3]))

##(ii) format to data.list
data.list<-fpca.format(data.m)
n<-length(data.list)    ##number of subjects
 if(n==0){ 
 print("error: no subject has more than one measurements!")
 return(0)
 }

##(iii)rescale time onto [0,1]
  data.list.new<-data.list
  for (i in 1:n){
  cur<-data.list[[i]][[1]]
  temp<-(cur[,2]-L1)/(L2-L1)
  temp[temp<0.00001]<-0.00001
  temp[temp>0.99999]<-0.99999 
  cur[,2]<-temp
  data.list.new[[i]][[1]]<-cur
  }
  
 ## (iv)estimate mean and substract
  #library(sm) 
  temp<-LocLin.mean(data.list.new,n,nmax, grids)  ##subtract estimated mean
  data.list.new<-temp[[1]]
  fitmu<-temp[[2]]

 ## (v)apply fpca.fit on all combinations of M and r in M.set and r.set and rescale back
  result<-NULL
  
  for (k in 1:length(r.set)){
  r.c<-r.set[k]
  print(paste("r=",r.c))
  result.c<-fpca.fit(M.set,r.c,data.list.new,n,nmax,grid.l,grids,iter.EM.num,sig.EM,ini.method,basis.method,
  sl.v,max.step,tol,cond.tol,M.self,no.EM,basis.EM)
                     
   if(is.vector(result.c[[1]])){
    print("warning: see warning code")
    } 

  ##rescale back
  grids.new<-grids*(L2-L1)+L1
  
  temp<-result.c$eigenvec
  M.set.u<-M.set[M.set>=r.c]   ##only when M>=r, there is a result
  if(length(M.set.u)==0){
  temp<-NULL
  }else{
  for(j in 1:length(M.set.u)){   
  temp[[j]]<-temp[[j]]/sqrt(L2-L1)
   }
  }
  
  result.c$eigenvec<-temp
  result.c$eigenval<-result.c$eigenval*(L2-L1)
  result.c$grids<-grids.new  

  result[[k]]<-result.c
  }  
  
  ##(vi) model selection
  mod.sele<-fpca.cv(result,M.set,r.set,tol=tol)
  temp<-mod.sele[[3]]
  index.r<-temp[1] 
  index.M<-temp[2]
  cv.result<-mod.sele[[1]]
  con.result<-mod.sele[[2]]
   
  ##(vii) return the selected model
  result.sele<-result[[index.r]]
  eigenf<-result.sele$eigenvec
  eigenf.sele<-eigenf[[index.M]][,,1]
  
  eigenv<-result.sele$eigenval
  eigenv.sele<-eigenv[1,,index.M]
  
  sig<-result.sele$sig
  sig.sele<-sig[1,index.M]
  
  temp.model<-c(M.set[index.M],r.set[index.r])
  names(temp.model)<-c("M","r")
  rownames(eigenf.sele)<-paste("eigenfunction",1:r.set[index.r])
  names(eigenv.sele)<-paste("eigenvalue",1:r.set[index.r])
  temp<-list("selected_model"=temp.model,"eigenfunctions"=eigenf.sele,"eigenvalues"=eigenv.sele,"error_var"=sig.sele^2,"fitted_mean"=fitmu,"grid"=grids.new,"cv_scores"=cv.result,"converge"=con.result)
 
 return(temp)  
}


                        






####################################################
############################ internal functions
##name:fpca.fit
##purpose: for a given dataset which is formatted in the form of data.list (e.g., by function fpca.format), 
##fit MLE by Newton for a set of  M (number of basis) and a fixed r (dimension of the process)

fpca.fit<-function(M.set,r,data.list,n,nmax,grid.l=seq(0,1,0.01),grids=seq(0,1,0.002),
iter.num=50,sig.EM=1,ini.method, basis.method="poly",sl.v,max.step=50,tol=1e-3,cond.tol=1e+10,
M.self=20,no.EM=FALSE,basis.EM="poly"){

##M.set--# of basis functions; r: dimension of the process
##data.list: data of the fpca.format format 
##n: sample size; nmax: maximum number of measurements per subject;
##grid.l: grids for loc; grids: denser grids for EM and Newton;
##iter.num: max number of iteration of EM; sig.EM: initial values of sig (erros sd.) for EM
##ini.method: initial method for Newton, one of "EM","loc","EM.self" (meaning using a fixed M.self); 
##basis.method: basis functions to use for Newton, one of "poly" (cubic Bspline), "ns" (nature spline)
##sl.v: shrinkage steps for Newton, e.g., sl.v<-c(0.5,0.5,0.5) means the first three steps are 0.5
##max.step: max number of iterations of Newton
##tol: tolerance level to determine convergence in Newton in terms of the l_2 norm of the gradient 
##cond.tol: tolerance to determine singularity in Newton in terms of condition number of the Hessian matrix
##M.self: dimension to use for EM.self; 
##no.EM: do EM (F) or not (T);
##basis.EM: basis to use for EM, one of "poly" or "ns"

 M.set<-M.set[M.set>=r]    ## only fit Newton for those M that is >= r.  
 M.l<-length(M.set)
  if(M.l==0){              ##if all M in the set M.set are less than r, return zero  
   print("all M<r")
   return(0)
  }

##results for return
 eigen.result<-array(0,dim=c(3,r,M.l))   ##estimated eigenvalues for three methods(Newton,ini, ini/EM) under different M
 eigenf.result<-NULL                     ##estimated eigenfunctions 
 sig.result<-matrix(0,3,M.l)             ##estimated error sd. for three methods under differen M
 like.result<-matrix(0,3,M.l)            ##likelihood for three methods under different M  
 cv.result <-numeric(M.l)                ##approximate CV score for Newton under different M
 converge.result<-numeric(M.l)           ##whether Newton converges under different M
 step.result<-numeric(M.l)               ## number of iterations for Newton to converge 
 cond.result<-numeric(M.l)               ## the condition number of the final Hessian in Newton
 
 if(!no.EM){                             ## if do EM 
 names.m<-c("newton","ini","EM")  
 }else{
 names.m<-c("newton","ini",ini.method)   ## if do not do EM
 }

 dimnames(eigen.result)<-list(names.m,NULL,NULL)
 rownames(sig.result)<-names.m
 rownames(like.result)<-names.m

###initial estimates by loc  
 if(ini.method=="loc"){ 
 IniVal<-try(Initial(r,ini.method,data.list,n,nmax,grid.l,grids))
   if(inherits(IniVal, "try-error")){
    print("warning: error in initial step by loc")
   return(-2)
   } 
   
 sig2hat<-IniVal[[1]]
 covmatrix.ini<-IniVal[[2]] 
 eigenf.ini<-IniVal[[3]]
 eigenv.ini<-IniVal[[4]]
 like.ini<-IniVal[[5]]
}

###initial values from EM with M=M.self
 if(ini.method=="EM.self"){ 
  IniVal<-try(Initial(r,ini.method="EM",data.list,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M.self,sig.EM=sig.EM))
  if(inherits(IniVal, "try-error")){
     print("warning: error in initial step by EM.self")
   return(-3)
   }

 sig2hat<-IniVal[[1]]
 covmatrix.ini<-IniVal[[2]] 
 eigenf.ini<-IniVal[[3]]
 eigenv.ini<-IniVal[[4]]
 like.ini<-IniVal[[5]]
}

 for(i in 1:M.l){
 #print(paste("M=", M.set[i]))
###################### EM
##(i)
 if(!no.EM){
 M.EM<-M.set[i]      ##number of basis use
 temp.EM<-try(EM(data.list,n,nmax,grids,M.EM,iter.num,r,basis.EM=basis.EM,sig.EM))
    if(inherits(temp.EM, "try-error")){
     print("warning: error in initial step by EM")
   return(-1)
   }

 EMeigenvec.est<-temp.EM[[1]]*sqrt(length(grids))
 EMeigenval.est<-temp.EM[[2]]
 EMsigma.est<-temp.EM[[3]]
 covmatrix.EM<-EMeigenvec.est%*%diag(EMeigenval.est)%*%t(EMeigenvec.est)

##(iii)likelihood of EM result
like.EM<-loglike.cov.all(covmatrix.EM,EMsigma.est,data.list,n)
 }else{ 
##use other initial value in the place of EM
 EMeigenvec.est<-eigenf.ini
 EMeigenval.est<-eigenv.ini
 EMsigma.est<-sqrt(sig2hat)
 like.EM<-like.ini 
 }



##########################newton's method
#####(0): Initial value for Newton's method

 if(ini.method=="EM"){
 sig2hat<-EMsigma.est^2
 covmatrix.ini<-covmatrix.EM
 eigenf.ini<-EMeigenvec.est
 eigenv.ini<-EMeigenval.est
 like.ini<-like.EM
}

####(i) projection on basis functions to get inital values of B and Lambda
##
M<-M.set[i]   ##number of basis to use
temp<-try(Proj.Ini(covmatrix.ini,M,r,basis.method,grids))
   if(inherits(temp, "try-error")){
     print("warning: error in Newton step due to projection")
    return(-4)
    }
B.ini<-temp[[1]]
Lambda.ini<-temp[[2]]        
 if(any(Lambda.ini<0)) {  ##should >0
   print("warning: error in Newton step: initial covariance not p.d.!")
   return(-4)
  }
  
if(sig2hat>0){
sig.ini<-sqrt(sig2hat)
}else{
sig.ini<-sqrt(Lambda.ini[r])
}

#### (ii)get auxilary things
if(basis.method=="poly"){ ##cubic B-spline with equal spaced knots
R.inv<-R.inverse(M,grids)
phi.aux<-apply(matrix(1:n),MARGIN=1,Phi.aux,M=M,basis.method=basis.method,data.list=data.list,R.inv=R.inv)
}

if(basis.method=="ns"){ ##nature spline
bs.u<-NS.orth(grids,df=M)/sqrt(grids[2]-grids[1])
phi.aux<-apply(matrix(1:n),MARGIN=1,Phi.aux,M=M,basis.method=basis.method,data.list=data.list,R.inv=bs.u, grid=grids)
}

#### (iv) Newton iterations
newton.result<-Newton.New(B.ini,phi.aux,sig.ini,Lambda.ini,data.list,n,sl.v,max.step,tol,cond.tol)
 step<-newton.result[[7]]
 like<-newton.result[[1]][1:step]
 B.up<-newton.result[[2]][,,step]
 Lambda.up<-newton.result[[3]][step,]
 sig.up<-newton.result[[4]][step]
 gradB<-newton.result[[5]][,,step]
 gradL<-newton.result[[6]][step,] 
 cond<-newton.result[[8]]
 error.c<-newton.result[[9]]

##check convergence
 converge<-max(abs(gradB),abs(gradL))          ##converged? if zero

#### (v) cross-validation score
if(error.c){
print("warning: cv not computed due to errors in the Newton iteration steps")
cv.score<-(-99)
}else{
cv.score<-try(CV(B.up,phi.aux,sig.up,Lambda.up,data.list,n))
error.c<-inherits(cv.score, "try-error")
if(error.c){
print("warning: cv not compuated due to errors in the CV step")
 cv.score<-(-99)
 }
}
   
   
   
## (vi)estimated eigenfunctions/values 
##eigenfunctions
eigenv.up<-Lambda.up[order(-Lambda.up)]
B.f<-B.up[,order(-Lambda.up)]
if(basis.method=="poly"){
eigenf.up<-t(apply(matrix(grids),1,Eigenf,basis.method=basis.method,B=B.f,R.inv=R.inv))
}
if(basis.method=="ns"){
eigenf.up<-t(apply(matrix(grids),1,Eigenf,basis.method=basis.method,B=B.f,R.inv=bs.u,grid=grids))
}

##
 dim.c<-c(r,length(grids),3)
 eigen.est<-array(0,dim.c)
 eigen.est[,,1]<-t(eigenf.up)        ##Newton
 eigen.est[,,2]<-t(eigenf.ini)       ##Initial
 eigen.est[,,3]<-t(EMeigenvec.est)   ##EM/initial   
 dimnames(eigen.est)<-list(NULL,NULL,names.m)
 eigenf.result[[i]]<-eigen.est
 
##eigenvalues, sig, like, cv, converge, step, cond
 eigen.temp<-rbind(eigenv.up,eigenv.ini,EMeigenval.est)  
 eigen.result[,,i]<-eigen.temp
 sig.result[,i]<-c(sig.up,sig.ini,EMsigma.est)
 like.result[,i]<-c(like[step],like.ini,like.EM)
 cv.result[i]<-cv.score
 converge.result[i]<-converge
 step.result[i]<-step
 cond.result[i]<-cond
}

####(vii) return result
result<-list("eigenval"=eigen.result,"sig"=sig.result,"-2*loglike"=like.result,
"cv"=cv.result,"converge"=converge.result,"step"=step.result,"condition#"=cond.result,
"eigenvec"=eigenf.result)
return(result)
}


##name: fpca.cv
##purpose: model selection based on the results of a set of models
 fpca.cv<-function(result.new,M.set,r.set,tol=1e-3){
  cv.mod<-matrix(-99,length(r.set),length(M.set))
  colnames(cv.mod)<-M.set
  rownames(cv.mod)<-r.set
 
  cond.mod<-cv.mod
  cv.sele<-cv.mod+1e+10

  for (k in 1:length(r.set)){
   r.c<-r.set[k]
   M.set.c<-M.set[M.set>=r.c]
   index.c<-sum(M.set<r.c)
   
    if(length(M.set.c)>0){
   for(j in 1:length(M.set.c)){
   cv.mod[k,j+index.c]<-result.new[[k]]$cv[j]
   cond.mod[k,j+index.c]<-result.new[[k]]$converge[j]
    if(cv.mod[k,j+index.c]!=(-99)&&cond.mod[k,j+index.c]<tol) 
     cv.sele[k,j+index.c]<-cv.mod[k,j+index.c] 
    }
   } 
   
  }
 
   index.r<-1
   index.M<-1
    
  for (j in 1:length(M.set)){
   for(k in 1:length(r.set)){
    if(cv.sele[k,j]<cv.sele[index.r,index.M]){
    index.r<-k
    index.M<-j
    }
   }
  } 
 
##
  temp<-c(index.r,index.M)
  names(temp)<-c("r","M")
  result<-list("cv"=cv.mod,"converge"=cond.mod,"selected model"=temp)   
  return(result)
 }





###name: fpca.format
##purpose: format the data into data.list as the input of fpca.fit and also exclude subjects with only one measurement 
fpca.format<-function(data.m){
##para: data.m: the data matrix with three columns: column 1: subject ID, column 2: observation, column 3: measurement time
##return: a list of n components where n is the number of subjects
##The ith component is a matrix of two columns: the first column is the measurements of the ith subject,
##and the second column is the corresponding times of measurements of the ith subject

ID<-unique(data.m[,1]) 
n<-length(ID)
data.list<-NULL

count<-1
  for(i in 1:n){
   N.c<-sum(data.m[,1]==ID[i]) 
   cur<-data.m[data.m[,1]==ID[i],]
   
   if(N.c>1){
   Obs.c<-as.numeric(cur[,2])
   T.c<-as.numeric(cur[,3])
   temp<-matrix(cbind(Obs.c,T.c),nrow=N.c, ncol=2)
   data.list[[count]]<-list(temp,NULL) 
   count<-count+1
   }
  
  }

 return(data.list)
 
 }



require("sm")
require("splines")


