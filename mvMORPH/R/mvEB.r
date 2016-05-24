################################################################################
##                                                                            ##
##                               mvMORPH: mvEB                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, geiger                                  ##
##                                                                            ##
################################################################################

mvEB<-function(tree, data, error=NULL, param=list(up=0), method=c("rpf","sparse","inverse","pseudoinverse","pic"), scale.height=FALSE, optimization=c("Nelder-Mead","L-BFGS-B","subplex"), control=list(maxit=20000), precalc=NULL, diagnostic=TRUE, echo=TRUE){

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
# select default model
method<-method[1]
# Check if there is missing cases
NA_val<-FALSE
Indice_NA<-NULL
if(any(is.na(data))){
    if(method!="pic" & method!="sparse"){
    NA_val<-TRUE
    }else{
     stop("NA values are allowed only with the \"rpf\",\"inverse\" or \"pseudoinverse\" methods")
    }
    Indice_NA<-which(is.na(as.vector(data)))
}
# bind error to a vector
if(!is.null(error)){error<-as.vector(error)}

# number of species (tip)
n<-dim(data)[1]
# number of variables
k<-dim(data)[2]
# method for the optimizer & algorithm
optimization<-optimization[1]
# scale height of the tree
if(scale.height==TRUE){
maxHeight<-max(nodeHeights(tree))
tree$edge.length<-tree$edge.length/maxHeight
}
# number of traits
k<-ncol(data)
if(is.null(k)){
k<-1
}

# number of parameters
npar<-(k*(k+1)/2)

if(is.null(param[["low"]])==TRUE){
    if(scale.height==FALSE) maxHeight<-max(nodeHeights(tree))
    low<-param$low<- log(10^-5)/maxHeight # Slater & Pennell 2013
   if(echo==TRUE) cat("No lower bound provided. Use of default setting \" ",low,"\"")

}else{low<-param$low}
if(is.null(param[["up"]])==TRUE){
  if(echo==TRUE)  cat("No upper bound provided. Use of default setting \"0\"")
    up<-param$up<-0
}else{up<-param$up}

if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
    
    tree<-precalc$tree
    C<-precalc$C1
    D<-precalc$D
    
    if(method=="sparse"){
        # Yale sparse format
        JAr<-precalc$JAr
        IAr<-precalc$IAr
        ch<-precalc$ch
        precalcMat<-precalc$V
    }
    
    
}else{
if(method!="pic"){
# compute the vcv
C<-vcv.phylo(tree) 
# Design matrix
D<-multD(tree,k,n,smean=TRUE)
    }else{
        C<-NULL
        D<-NULL
    }

if(method=="sparse"){
    V<-kronecker((matrix(1,k,k)+diag(k)), C)
    # spam object
    precalcMat<-as.spam(V);
    # precal the cholesky
    if(is.null(param[["pivot"]])){pivot<-"MMD"}else{pivot<-param$pivot}
    ch<-chol(precalcMat,pivot=pivot)
    # Yale Sparse Format indices
    JAr<-precalcMat@colindices-1
    IAr<-precalcMat@rowpointers-1
}else{
    ch<-NULL
    precalcMat<-NULL
}

}
##--------------Likelihood_functions--------------------------------------------------##


# switch method
switch(method,
"rpf"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc){
        V<-.Call("kroneckerEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n))
        loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA) ######### A modif pour sparse
        return(loglik)
    }
},
"sparse"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc){
        V<-.Call("kroneckerSparEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n),  IA=as.integer(precalcMat@rowpointers-1), JA=as.integer(precalcMat@colindices-1), A=precalcMat@entries)
        loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=NULL) ######### A modif pour sparse
        return(loglik)
    }
},
"pseudoinverse"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc){
        V<-.Call("kroneckerEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n))
        loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA) ######### A modif pour sparse
        return(loglik)
    }
},
"inverse"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc){
    V<-.Call("kroneckerEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA) ######### A modif pour sparse
    return(loglik)
    }
},
"pic"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc){
        res<-.Call("PIC_gen", x=dat, n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=list(tree$edge.length), times=tt, rate=as.double(rep(beta,k)), Tmax=1, Model=as.integer(1), mu=1, sigma=sig)
        logl<- -0.5 * ( n * k * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
       return(list(logl=logl,anc=res[[7]], sigma=res[[2]]))
    }
}
)



##--------------Maximum Likelihood estimation of EB---------------------------##  

likEB <- function(dat, error, C, D, beta, sig, k, n, method, precalc, precalcMat) {

    loglik<-eb_fun_matrix(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc) ######### A modif pour sparse

    list(loglik=-loglik$logl, ancstate=loglik$anc)
}


##---------------Optimizing function------------------------------------------##
# initial value for exponent parameter
if(is.null(param[["beta"]])==TRUE){
ebval<-runif(1,param$low,param$up)*0.1
}else{
ebval<-param$beta
}

# initial values for the optimizer
if(is.null(param[["sigma"]])==TRUE){
sig1<-varBM(tree,data,n,k)
sig1<-sym.unpar(sig1)
}else{
sig1<-param$sigma
}

if(method=="pic"){
    warning<-eval_polytom(tree)
    tree<-reorder.phylo(tree,"postorder")
    # times from the root
    tt<-.Call("times_root", brlength=tree$edge.length, edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), ntip=as.integer(n), Nnode=as.integer(tree$Nnode))
    
    if(optimization=="subplex"){
        estim <- subplex(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sym.par(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, hessian=TRUE, control=control)
    }else{
        estim <- optim(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sym.par(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, gr=NULL, hessian=TRUE, method = optimization, control=control)
    }
}else{


if(optimization=="subplex"){
    estim <- subplex(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sym.par(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, hessian=TRUE, control=control)
   }else{
    estim <- optim(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sym.par(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, gr=NULL, hessian=TRUE, method = optimization, control=control)
   }
}
##-----------------Summarizing results----------------------------------------##
# sigma matrix
resultList<-sym.par(estim$par[1+seq_len(npar)])
colnames(resultList)<-colnames(data)
rownames(resultList)<-colnames(data)
#ancestral states estimates
anc<-matrix(likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,estim$par[1]), sig=sym.par(estim$par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$ancstate,nrow=1)
rownames(anc)<-"theta"
colnames(anc)<-colnames(data)

# rate parameter
r=ratevalue(up,low,estim$par[1])
# LogLikelihood
LL<--estim$value
# models parameters
nparam=k+npar+1 # length(estim$par) # AIC
AIC<--2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1989
##---------------------Diagnostics--------------------------------------------##

if(estim$convergence==0 & diagnostic==TRUE){  
cat("\n","successful convergence of the optimizer","\n") 
}else if(estim$convergence==1 & diagnostic==TRUE){  
cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
}else if(diagnostic==TRUE){  
cat("\n","convergence of the optimizer has not been reached, try simpler model","\n") 
}
# Hessian eigen decomposition to check the derivatives
hess<-eigen(estim$hessian)$values
if(any(hess<0)){
hess.value<-1
if(diagnostic==TRUE){
cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
}else{
hess.value<-0
if(diagnostic==TRUE){
cat("a reliable solution has been reached","\n")}
}

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("-- Summary results for Early Burst or ACDC model --","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters","\n")
cat("Rate change:","\n")
cat("______________________","\n")
print(r)
cat("\n")
cat("Estimated rates matrix","\n")
cat("______________________","\n")
print(resultList)
cat("\n")
cat("Estimated root states","\n")
cat("______________________","\n")
print(anc)
cat("\n")
}
##-------------------Save infos in parameters---------------------------------##
param$nparam<-nparam
param$nbspecies<-n
param$ntraits<-k
param$nregimes<-1
param$method<-method
param$optimization<-optimization
param$traits<-colnames(data)
param$model<-if(r>0){"AC"}else if(r<0){"EB"}else{"ACDC"}
##-------------------Store results--------------------------------------------##
 
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=anc, beta=r, sigma=resultList,  convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.value, param=param)
class(results)<-c("mvmorph","acdc")
invisible(results)

}
