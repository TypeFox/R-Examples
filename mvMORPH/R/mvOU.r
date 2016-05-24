################################################################################
##                                                                            ##
##                               mvMORPH: mvOU                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex                                 ##
##                                                                            ##
################################################################################

mvOU<-function(tree,data,error=NULL,model=c("OUM","OU1"),param=list(sigma=NULL,alpha=NULL, vcv="mvmorph", decomp=c("symmetricPositive","symmetric","nsymPositive","nsymmetric","diagonal","upper","lower")),method=c("rpf","sparse","inverse","pseudoinverse","univarpf"),scale.height=FALSE, optimization=c("L-BFGS-B","Nelder-Mead","subplex"),control=list(maxit=20000),precalc=NULL,diagnostic=TRUE, echo=TRUE){

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
# choose method for the likelihood computation
method=method[1]
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

# Check the order of the dataset and the phylogeny 

	if(!is.null(rownames(data))) { 
        if(any(tree$tip.label==rownames(data))){
	data<-data[tree$tip.label,] }else if(echo==TRUE){
    cat("row names of the data matrix must match tip names of your phylogeny!","\n")
	}}else if(echo==TRUE){
	cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
	}
# Check the tree
    if(is.null(tree[["mapped.edge"]])==TRUE){
        model<-"OU1"
      if(echo==TRUE)  cat("No selective regimes mapped on the tree, only a OU1 model could be estimated","\n")
    }

##-------------------------Calculation of parameters--------------------------##
   
  # number of species
n<-length(tree$tip.label)
  # choose model
model=model[1]

  # choose the optimization method
optimization=optimization[1]
  # Pull (alpha) matrix decomposition
if(is.null(param[["decomp"]])){
    decomp<-param$decomp<-"symmetric" ## for "ouch" vcv select "symmetricPositive" instead
}else{
    decomp<-param$decomp[1]
}
  # option for computing the variance covariance matrix
  if(is.null(param[["vcv"]])==TRUE){
      if(method=="sparse"){
          vcvtype<-"sparse"
      }else{
      if(is.ultrametric(tree)==TRUE & decomp=="symmetricPositive"){
        vcvtype<-"ouch"
      if(echo==TRUE)   cat("The tree is ultrametric, the \"ouch\" VCV is used as default setting with symmetric positive definite alpha matrix. See ?mvOU","\n")
      }else{
        vcvtype<-"mvmorph"
       if(echo==TRUE) cat("The tree is not ultrametric or the alpha matrix is not constrained to be symmetric positive definite, the \"mvmorph\" VCV is thus used. See ?mvOU","\n")
      }
    }
      
  }else{
      vcvtype<-param$vcv
  }
  
  if(vcvtype!="sparse" & method=="sparse"){
      vcvtype<-"sparse"
   if(echo==TRUE)   cat("Only \"sparse\" VCV could be used with the \"sparse\" method. See ?mvOU","\n")
      #method<-"sparse"
  }
  if(vcvtype=="ouch" & decomp!="symmetricPositive"){
      decomp<-param$decomp<-"symmetricPositive"
     if(echo==TRUE) cat("Only \"symmetricPositive\" alpha matrices could be used with the \"ouch\" vcv method. See ?mvOU","\n")
      # otherwise the program crash
  }
  
  # root estimation
  if(is.null(param[["root"]])!=TRUE){
      if(param[["root"]]==TRUE || param[["root"]]==FALSE || param[["root"]]=="stationary"){
     root<-param$root
        }else{
      stop("Only TRUE,FALSE or \"stationary\" are accepted for the \"root\" argument in \"param\"")
        }
        if(param[["root"]]=="stationary" & model=="OU1"){
     root<-param$root<-FALSE
        }
  }else if(vcvtype=="ouch" | vcvtype=="sparse"){
      root<-FALSE
  }else if(vcvtype=="mvmorph" & is.ultrametric(tree)){
      root<-FALSE
  }else{
      root<-TRUE
  }
  

  # regimes number
  if(model=="OUM"){
      k<-length(colnames(tree$mapped.edge))
  }else{ k<-1 }
  # number of traits
  if(is.vector(data)){
  p<-1 }else{ p<-ncol(data)}
  
  # if univariate
  if(p==1){
      
if(is.null(param[["alpha"]])!=TRUE){
    if(length(param$alpha)>1){
    stop("You must provide a unique starting value for univariate models")
    }
    if(!is.numeric(param$alpha[[1]])){
    stop("You can't use a constraint on univariate models")
    }
alpha<-param$alpha
}

if(is.null(param[["sigma"]])!=TRUE){
    if(length(param$sigma)>1){
    stop("You must provide a unique starting value for univariate models")
    }
    if(!is.numeric(param$sigma[[1]])){
    stop("You can't use a constraint on univariate models")
    }
sigma<-param$sigma
}

if(vcvtype=="mvmorph"){ # a modifier pour integrer sparse et rpf?
    # si explicite ->
    if(method=="univarpf"){vcvtype<-"univarpf"}else{ vcvtype<-"univar"}
    }
decomp<-param$decomp<-"symmetricPositive"
}
  if(method=="univarpf"){
      if(p!=1){ method<-"rpf" }
      if(vcvtype=="ouch"){ vcvtype<-"univarpf"
        if(echo==TRUE)    cat("Only \"mvmorph\" VCV could be used with the \"univarpf\" method. See ?mvOU","\n")
      }
  }

  # bind data to a vector
if(is.matrix(data)){
dat<-as.vector(data) }else{ dat<-as.vector(as.matrix(data))}
  # bind error to a vector
if(!is.null(error)){error<-as.vector(error)}
## initial alpha and sigma matrix if not provided / or constrained models
# sigma matrix  
constrSigma<-FALSE 
if(is.null(param[["sigma"]])){
 sigma<-param$sigma<-sym.unpar(diag(0.1,p))
}
if(!is.numeric(param$sigma[[1]])){
    if(length(param$sigma)==2){
    if(is.numeric(param$sigma[[2]])){sigma=param$sigma[[2]]}
    }else{ sigma=rep(0.1,p) }
  constrSigma<-TRUE
}else{
    sigma<-param$sigma
}


# alpha matrix
if(is.null(param[["alpha"]])){
    if(decomp=="symmetric" | decomp=="symmetricPositive" | decomp=="upper" | decomp=="lower"){
        alpha<-param$alpha<-sym.unpar(diag(runif(p),p))
    }else if(decomp=="nsymmetric" | decomp=="nsymPositive"){
        alpha<-param$alpha<-runif(p*p)
    }else if(decomp=="diagonal" | decomp=="diagonalPositive" | decomp=="univariate"){
        alpha<-param$alpha<-runif(p)
    }
}
if(!is.numeric(param$alpha[[1]])){
if(length(param$alpha)==2){
  if(is.numeric(param$alpha[[2]])){
  alpha=param$alpha[[2]]}}else{ alpha=rep(0.1,p)}
   decomp<-"diagonal"
}else{
    alpha<-param$alpha
}
     
# number of parameters for each matrix
nalpha<-length(alpha)
nsigma<-length(sigma)

# function for alpha parameterization
buildA<-function(x,p){matrixParam(x,p,decomp)}

# number of columns of the design matrix
if(root==TRUE){sizeD<-(k+1)*p}else{sizeD<-k*p}

# Define the variance-covariance and weights matrix for the likelihood function

switch(vcvtype,
"ouch"={
    ou_fun_matrix<-function(bt,n,alphaA,sigmA,epochs,listReg,mod_stand){
        V<-.Call("simmap_covar", as.integer(n), bt=bt, lambda=alphaA$values, S=alphaA$vectors, sigmasq=sigmA)
        W<-.Call("simmap_weights",nterm=as.integer(n), epochs=epochs,lambda=alphaA$values,S=alphaA$vectors,beta=listReg)
        list(V=V, W=W)
    }
},
"sparse"={
     ou_fun_matrix<-function(bt,n,alphaA,sigmA,epochs,listReg,mod_stand){
         V<-.Call("mvmorph_covar_ou_sparse", A=as.double(precalcMat@entries), JA=as.integer(JAr), IA=as.integer(IAr), as.integer(n), bt=bt, lambda=alphaA$values, S=alphaA$vectors, sigmasq=sigmA, S1=alphaA$invectors)
         W<-.Call("mvmorph_weights",nterm=as.integer(n), epochs=epochs,lambda=alphaA$values,S=alphaA$vectors,S1=alphaA$invectors,beta=listReg,root=as.integer(mod_stand))
         list(V=V, W=W)
     }
},
"mvmorph"={
     ou_fun_matrix<-function(bt,n,alphaA,sigmA,epochs,listReg,mod_stand){
         V<-.Call("mvmorph_covar_mat", as.integer(n), bt=bt, lambda=alphaA$values, S=alphaA$vectors, sigmasq=sigmA, S1=alphaA$invectors)
         W<-.Call("mvmorph_weights",nterm=as.integer(n), epochs=epochs,lambda=alphaA$values,S=alphaA$vectors,S1=alphaA$invectors,beta=listReg,root=as.integer(mod_stand))
         list(V=V, W=W)
     }
},
"univar"={
     ou_fun_matrix<-function(bt,n,alphaA,sigmA,epochs,listReg,mod_stand){
         V<-.Call("mvmorph_covar_ou",A=bt,alpha=alphaA$values, sigma=sigmA)
         W<-.Call("mvmorph_weights",nterm=as.integer(n), epochs=epochs,lambda=alphaA$values,S=1,S1=1,beta=listReg,root=as.integer(mod_stand))
         list(V=V, W=W)
     }
},
"univarpf"={
     ou_fun_matrix<-function(bt,n,alphaA,sigmA,epochs,listReg,mod_stand){
         V<-.Call("mvmorph_covar_ou_rpf",A=mt$mDist,alpha=alphaA$values, sigma=sigmA)
         W<-.Call("mvmorph_weights",nterm=as.integer(n), epochs=epochs,lambda=alphaA$values,S=1,S1=1,beta=listReg, root=as.integer(mod_stand))
    # time gain only for very huge phylogeny
    list(V=V, W=W)
     }
})


##-----------------------Precalc-on-------------------------------------------##
# precalc
if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
    tree<-precalc$tree
    mt<-list(mDist=precalc$C1)
    root<-precalc$param$root
    if(method=="sparse"){
        # Yale sparse format
        JAr<-precalc$JAr
        IAr<-precalc$IAr
        ch<-precalc$ch
        precalcMat<-precalc$V
    }
    
    listReg<-precalc$listReg
    epochs<-precalc$epochs
    model<-precalc$model
    if(model=="OU1"){
        k<-1
    }else{
        k<-length(colnames(tree$mapped.edge))
    }
    
    # number of columns of the design matrix
    if(root==TRUE){sizeD<-(k+1)*p}else{sizeD<-k*p}
    if(root==FALSE){
        mod_stand<-0 # the root is not provided nor assumed to be one of the selected regimes, so we rowstandardize (could be optional)
    }else if(root==TRUE){
        k<-k+1
        mod_stand<-0
    }else if(root=="stationary"){
        mod_stand<-1
    }
    
}else{
    
# number of columns of the design matrix
if(root==TRUE){sizeD<-(k+1)*p}else{sizeD<-k*p}
# ancestor times calculation
mt<-mTime(tree,scale.height)
# max node height for standardisation
mb<-max(nodeHeights(tree))

    if(method=="sparse"){
        V<-kronecker((matrix(1,p,p)+diag(p)), mt$mDist)
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


##-----------------------Precalculate regime indexation-----------------------##
## a mettre dans pre-calc aussi?

  # root to tip lineage indexation
root2tip <- .Call("seq_root2tipM", tree$edge, n, tree$Nnode)
# Si OU1 sur un objet 'phylo'
if(model=="OU1"){
if(scale.height==TRUE){
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$edge.length[val]<-tree$edge.length[val]/mb},simplify=FALSE)))
} ,simplify=FALSE)}else{
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$edge.length[val]<-tree$edge.length[val]},simplify=FALSE)))
} ,simplify=FALSE)
}
}else{
# Donnees de temps par regimes et par branches
if(scale.height==TRUE){
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$maps[[val]]<-tree$maps[[val]]/mb},simplify=FALSE)))
} ,simplify=FALSE)
}else{
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$maps[[val]]<-tree$maps[[val]]},simplify=FALSE)))
} ,simplify=FALSE)
}
}

# Indexer les regimes
if(model=="OUM"){ 
# Indexing factors
if(root==FALSE | root=="stationary"){
    facInd<-factor(colnames(tree$mapped.edge))
    }else if(root==TRUE){
        facInd<-factor(c("_root_state",colnames(tree$mapped.edge)))
    }

    indice<-lapply(1:n,function(z){rev(unlist(lapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); factor(names(tree$maps[[val]]),levels=facInd)})))})

}else if(model=="OU1"){
    
    if(root==TRUE){
    facInd<-factor(c("_root_state","theta_1"))
    indice<-lapply(1:n,function(z){ as.factor(rep(facInd[facInd=="theta_1"],length(valLineage[[z]])))})
    }else{
    indice<-lapply(1:n,function(z){ as.factor(rep(1,length(valLineage[[z]])))})
    }
}
# Liste avec dummy matrix
if(root==FALSE){
    indiceA<-indiceReg(n,indice, facInd, FALSE)
    mod_stand<-0 # the root is not provided nor assumed to be one of the selected regimes, so we rowstandardize (could be optional)
}else if(root==TRUE){
    indiceA<-indiceReg(n,indice, facInd, TRUE)
    k<-k+1
    mod_stand<-0
}else if(root=="stationary"){
    indiceA<-indiceReg(n,indice, facInd, FALSE)
    mod_stand<-1
}


listReg<-sapply(1:n,function(x){sapply(1:p,function(db){regimeList(indiceA[[x]],k=k,root)},simplify=FALSE)},simplify=FALSE)


# mapped epochs
epochs<-sapply(1:n,function(x){lineage<-as.numeric(c(cumsum(valLineage[[x]])[length(valLineage[[x]])],(cumsum(valLineage[[x]])[length(valLineage[[x]])]-cumsum(valLineage[[x]])))); lineage[which(abs(lineage)<1e-15)]<-0; lineage },simplify=FALSE)

}# end of precalc option

##-----------------------Likelihood Calculation-------------------------------##

devianc<-function(alpha,sigma,dat,error,mt){
    
     alphaA<-buildA(alpha,p)
     if(constrSigma==TRUE){sigmA<-diag(diag(sigma%*%t(sigma)))}else{sigmA<-sym.par(sigma)}  #  trick... il ne faut pas de valeurs negatives sous risque de faire crasher
     
     matEstim<-ou_fun_matrix(bt=mt$mDist,n=n,alphaA,sigmA,epochs,listReg,mod_stand)

     # if (any(is.nan(diag(matEstim$V))) || any(is.infinite(diag(matEstim$V)))) return(1000000)
  
     loglik<-loglik_mvmorph(dat,matEstim$V,matEstim$W,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)
    

		if(is.infinite(loglik$logl)){
			return(10000000)
		}
        
		return(-loglik$logl)
	
	}
##----------------------Likelihood optimization-------------------------------##
if(optimization!="subplex"){
estim <- optim(par=c(alpha,sigma),fn = function (par) { devianc(alpha=par[seq_len(nalpha)],sigma=par[nalpha+seq_len(nsigma)],error=error,dat=dat,mt=mt)},gr=NULL,hessian=TRUE,method=optimization,control=control)
}else{
estim <- subplex(par=c(alpha,sigma),fn = function (par) {devianc(alpha=par[seq_len(nalpha)],sigma=par[nalpha+seq_len(nsigma)],error=error,dat=dat,mt=mt)},hessian=TRUE,control=control)
}
                   
##---------------------theta estimation---------------------------------------##

est.theta<-function(estimML){

    alphaA<-buildA(estimML[seq_len(nalpha)],p)
    if(constrSigma==TRUE){sigmA<-diag(diag(tcrossprod(estimML[nalpha+seq_len(nsigma)])))}else{sigmA<-sym.par(estimML[nalpha+seq_len(nsigma)])}

    matEstim<-ou_fun_matrix(bt=mt$mDist,n=n,alphaA,sigmA,epochs,listReg,mod_stand)
    V<-matEstim$V
    W<-matEstim$W

    mvmorphEstim<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method,ch=ch,precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)
		
    list(theta=mvmorphEstim$anc, V=V, W=W)
	}

res.theta<-est.theta(estim$par)$theta

##---------------------Diagnostics--------------------------------------------##
hess<-eigen(estim$hessian)$values

if(estim$convergence==0 & diagnostic==TRUE){ 
cat("successful convergence of the optimizer","\n")
}else if(estim$convergence==1 & diagnostic==TRUE){  
cat("maximum limit iteration has been reached, please consider increase maxit","\n")
}else if(diagnostic==TRUE){
cat("convergence of the optimizer has not been reached, try simpler model","\n")
}

if(any(hess<0)){
hess.val<-1  
if(diagnostic==TRUE){
cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
}else{
hess.val<-0
if(diagnostic==TRUE){
cat("a reliable solution has been reached","\n")}
}

##-------------------Summarize Results----------------------------------------##
LL<- -estim$value 
nparam=nalpha+nsigma+(k*p)

# maximum likelihood estimates of alpha and sigma
estim.alpha<-estim$par[seq_len(nalpha)]
estim.sigma<-estim$par[nalpha+seq_len(nsigma)]

# alpha matrix
alpha.mat<-buildA(estim.alpha,p)$A
# sigma matrix
if(constrSigma==TRUE){
sigma.mat<-diag(diag(estim.sigma%*%t(estim.sigma)))
}else{sigma.mat<-sym.par(estim.sigma)}
colnames(sigma.mat)<-rownames(sigma.mat)<-colnames(data)
colnames(alpha.mat)<-rownames(alpha.mat)<-colnames(data)
# AIC
AIC<- -2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) # Hurvich et Tsai, 1989
# matrix of estimated theta values
theta.mat<-matrix(res.theta,k)
if(model=="OUM"){
    if(root==TRUE){
     rownames(theta.mat)<-c("theta_0",colnames(tree$mapped.edge))
    }else{
     rownames(theta.mat)<-colnames(tree$mapped.edge)
    }
}else{
    if(root==TRUE){
    rownames(theta.mat)<-c("theta_0","theta_1")
    }else{
    rownames(theta.mat)<-"OU1"
    }
}
colnames(theta.mat)<-colnames(data)

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("-- Summary results --","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters","\n")
cat("\n")
cat("Estimated theta values","\n")
cat("______________________","\n")
print(theta.mat)
cat("\n")
cat("ML alpha values","\n")
cat("______________________","\n")
print(alpha.mat)
cat("\n")
cat("ML sigma values","\n")
cat("______________________","\n")
print(sigma.mat)
}

##-------------------Save infos in parameters---------------------------------##
param$nparam<-nparam
param$nbspecies<-n
param$ntraits<-p
param$nregimes<-k
param$method<-method
param$optimization<-optimization
param$traits<-colnames(data)
param$names_regimes<-colnames(tree$mapped.edge)
param$model<-model
param$root<-root
param$vcv<-vcvtype
param$decomp<-decomp

##------------------List results----------------------------------------------##


results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha=alpha.mat, sigma=sigma.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
class(results)<-c("mvmorph","ou")
invisible(results)

}


halflife <- function(object) UseMethod("halflife")
stationary <- function(object) UseMethod("stationary")
