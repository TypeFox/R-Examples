################################################################################
##                                                                            ##
##                               mvMORPH: mvSHIFT                             ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, subplex                                 ##
##                                                                            ##
################################################################################

mvSHIFT<-function(tree,data,error=NULL,param=list(age=NULL,sigma=NULL,alpha=NULL,sig=NULL,beta=NULL),model=c("ER","RR","EC","RC","SR","EBOU","OUEB","EBBM","BMEB"), method=c("rpf","sparse","inverse","pseudoinverse"),scale.height=FALSE, optimization=c("L-BFGS-B","Nelder-Mead","subplex"),control=list(maxit=20000),precalc=NULL,diagnostic=TRUE,echo=TRUE){
   
    #set data as a matrix if a vector is provided instead
    if(!is.matrix(data)){data<-as.matrix(data)}
    # choose method for the likelihood computation
    method=method[1]
    if(method!="rpf" & method!="inverse" & method!="pseudoinverse" & method!="sparse"){
        stop("Please choose either the \"rpf\", \"inverse\", \"pseudoinverse\", or \"sparse\" method")
    }
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
    # number of species (tip)
    n<-dim(data)[1]
    # number of variables
    p<-dim(data)[2]
    # choose model
    model=model[1]
    # choose a method for the optimizer
    optimization=optimization[1]
    if(is.null(param[["up"]])){up<-param$up<-2}else{up<-param$up}
    if(is.null(param[["low"]])){low<-param$low<- -2}else{low<-param$low}
    #set age shift with make.era from phytools if the age is provided
    if(!is.null(param[["age"]])){
    tot<-max(nodeHeights(tree))
    limit=tot-param$age
    tree<-make.era.map(tree, c(0,limit))
    }
    #Shift rate model
    if(model=="SR"){
        mvBM(tree,data,model="BMM",param=param,scale.height=scale.height,diagnostic=diagnostic,method=method,optimization=optimization,echo=echo,control=control,error=error)
    }else{
##---------------------------Precalc-parameters-----------------------------##

if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
    tree<-precalc$tree
    if(is.null(tree[["mapped.edge"]])){
    stop("The tree is not in SIMMAP format with mapped time slices")
    }
    if(!is.null(param[["age"]])){
  if(echo==TRUE)  cat("The tree specified in precalc is assumed to represent a shift at :",param$age)
    }
    # multiple vcv list
    vcvList<-precalc$C2
    # weight matrix
    W<-precalc$D
    # Yale sparse format
    JAr<-precalc$JAr
    IAr<-precalc$IAr
    ch<-precalc$ch
    precalcMat<-precalc$V
    
}else{
    # scale tree
    if(scale.height==TRUE){
        maxHeight<-max(nodeHeights(tree))
        tree$edge.length<-tree$edge.length/maxHeight
        tree$mapped.edge<-tree$mapped.edge/maxHeight
    }
    # compute the design matrix
    W <-multD(tree,p,n,smean=TRUE)
  
    # compute the variance-covariances matrices for each map
    multi.tre<-list()
    class(multi.tre)<-"multiPhylo"
	vcvList<-list()
	for(i in 1:ncol(tree$mapped.edge)){
		multi.tre[[i]]<-tree
		multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		if(any(tree$tip.label==rownames(data))) { 
		vcvList[[i]]<-temp[rownames(data),rownames(data)]
		}else{
		vcvList[[i]]<-temp
		}
	}
    if(method=="sparse"){
        ## sparse precalcs
        #require(spam)
        # Temporary multivariate VCV
        V<-kronecker((matrix(1,p,p)+diag(p)), vcv.phylo(tree))
        # spam object
        precalcMat<-as.spam(V);
        # precal the cholesky
        if(is.null(param[["pivot"]])){pivot<-"MMD"}else{pivot<-param$pivot}
        ch<-chol(precalcMat,pivot=pivot)
        # Yale Sparse Format indices
        JAr<-precalcMat@colindices-1
        IAr<-precalcMat@rowpointers-1
        # list precalc
    }
}

##---------------------------Define initial values--------------------------##
# Pull (alpha) matrix decomposition
if(is.null(param[["decomp"]])){
    decomp<-param$decomp<-"symmetricPositive"
}else{
    decomp<-param$decomp[1]
}
# initial sigma matrix if not provided
if(is.null(param[["sigma"]])){sigma=sym.unpar(diag(1,p))}else{sigma=param$sigma}
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
# second sigma matrix
if(is.null(param[["sig"]])){sig=sym.unpar(diag(1,p))}else{sig=param$sig}
# beta matrix
if(is.null(param[["beta"]])){beta=0.01}else{beta=param$beta} # 1 beta cf. mvEB

nalpha<-length(alpha)
nsigma<-length(sigma)

    #set parameters for ecological constraint model:
    if(model=="EC"|| model=="RC" || model=="BMOU" || model=="BMOUi"){
        before<-2
        after<-1
        if(model=="EC"||model=="BMOU"){
            model<-"ER"
            nmod<-"Ecological Constraint model"
            nsig<-nalpha
        }else if(model=="RC" || model=="BMOUi"){
            model<-"RR"
            nmod<-"Radiation and Ecological Constraint model"
            nsig<-nalpha+nsigma
        }
    }else if(model=="ER"|| model=="RR" || model=="OUBM" || model=="OUBMi"){
        before<-1
        after<-2
        if(model=="RR" || model=="OUBMi"){
            model<-"RR"
            nmod<-"Release and Radiate model"
            nsig<-nalpha+nsigma
        }else if(model=="ER" || model=="OUBM"){
            model<-"ER"
            nmod<-"Ecological Release model"
            nsig<-nalpha
        }
    }else if(model=="EBOU" || model=="EBOUi"){# OU-ACDC models
        before<-2
        after<-1
        if(model=="EBOU"){
            model<-"OV"
            nmod<-"ACDC to OU process with the same rates"
            nsig<-nalpha
            nbeta<-nsig+nsigma
        }else if(model=="EBOUi"){
            model<-"OVG"
            nmod<-"ACDC to OU process with independent rates"
            nsig<-nalpha+nsigma
            nbeta<-nsig+nsigma
        }
    }else if(model=="OUEB" || model=="OUEBi"){
        before<-1
        after<-2
        if(model=="OUEB"){
            model<-"OV"
            nmod<-"OU to ACDC process with the same rates"
            nsig<-nalpha
            nbeta<-nsig+nsigma
        }else if(model=="OUEBi"){
            model<-"OVG"
            nmod<-"OU to ACDC process with independent rates"
            nsig<-nalpha+nsigma
            nbeta<-nsig+nsigma
        }
    }else if(model=="EBBM" || model=="EBBMi"){ # BM-ACDC models
        before<-2
        after<-1
        if(model=="EBBM"){
            model<-"CV"
            nmod<-"ACDC to BM process with the same rates"
            nalpha<-1
            nbeta<-1
            nsig<-nbeta
        }else if(model=="EBBMi"){
            model<-"CVG"
            nmod<-"ACDC to BM process with independent rates"
            nalpha<-1
            nbeta<-1
            nsig<-nbeta+nsigma
        }
    }else if(model=="BMEB" || model=="BMEBi"){
        before<-1
        after<-2
        if(model=="BMEB"){
            model<-"CV"
            nmod<-"BM to ACDC process with the same rates"
            nalpha<-1
            nbeta<-1
            nsig<-nbeta
        }else if(model=="BMEBi"){
            model<-"CVG"
            nmod<-"BM to ACDC process with independent rates"
            nalpha<-1
            nbeta<-1
            nsig<-nbeta+nsigma
        }
    }


  # function for alpha parameterization
  buildA<-function(x,p){matrixParam(x,p,decomp)}

  # bind error to a vector
  if(!is.null(error)){error<-as.vector(error)}
  
  # bind data to a vector
  if(is.matrix(data)){ dat<-as.vector(data) }else{ dat<-as.vector(as.matrix(data))}
  
  # number of columns of the design matrix
  sizeD<-ncol(W)
##------------------------Likelihood function---------------------------------##
# ou-bm
lik.shift_oubm<-function(alpha,sigma,sig,dat,error,vcvList,p){
  if(p==1){
  Vou<-.Call("mvmorph_covar_ou",A=vcvList[[before]],alpha=alpha$values, sigma=sigma)
  }else{
  Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=vcvList[[before]],lambda=alpha$values,S=alpha$vectors,sigmasq=sigma, S1=alpha$invectors)
  }
  V<-.Call("kronecker_shift", R=sig, C=vcvList[[after]], Rrows=as.integer(p), Crows=as.integer(n), V=Vou)
  
  loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)

     list(loglik=-loglik$logl, theta=loglik$anc)
}

# eb-ou
lik.shift_ebou<-function(alpha,sigma,sig,beta,dat,error,vcvList,p){
    if(p==1){
        Vou<-.Call("mvmorph_covar_ou",A=vcvList[[before]],alpha=alpha$values, sigma=sigma)
    }else{
        Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=vcvList[[before]],lambda=alpha$values,S=alpha$vectors,sigmasq=sigma, S1=alpha$invectors)
    }
    V<-.Call("kronecker_shiftEB_OU", R=sig, C=vcvList[[after]], beta=matrix(beta,p,p), Rrows=as.integer(p),  Crows=as.integer(n), V=Vou)
    
    loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)
    
       list(loglik=-loglik$logl, theta=loglik$anc)
}

# eb-bm
lik.shift_ebbm<-function(beta,sigma,sig,dat,error,vcvList,p){
    # C2 is the EB matrix
    V<-.Call("kronecker_shiftEB_BM", R1=sig, R2=sigma, C1=vcvList[[before]], C2=vcvList[[after]], beta=matrix(beta,p,p), Rrows=as.integer(p),  Crows=as.integer(n))
    
    loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)
    
    list(loglik=-loglik$logl, theta=loglik$anc)
}

# sparse matrix ou-bm
lik.shift_oubm_sparse<-function(alpha,sigma,sig,dat,error,vcvList,p){
    Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=vcvList[[before]],lambda=alpha$values,S=alpha$vectors,sigmasq=sigma, S1=alpha$invectors)
    V<-.Call("kroneckerSpar_shift", R=sig, C=vcvList[[after]], Rrows=as.integer(p),  Crows=as.integer(n), V=Vou, IA=as.integer(IAr), JA=as.integer(JAr), A=as.double(precalcMat@entries))
    
    loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method="sparse",ch=ch,precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)
    
       list(loglik=-loglik$logl, theta=loglik$anc)
}

# sparse matrix eb-ou
#lik.shift_ebou_sparse<-function(alpha,sigma,sig,beta,dat,error,vcvList,p){
#    eig<-eigen(alpha)
#    Vou<-.Call("mvmorph_covar_ou_sparse", A=as.double(precalcMat@entries), JA=as.integer(JAr), IA=as.integer(IAr), as.integer(n*p), bt=vcvList[[before]], lambda=eig$values, S=eig$vectors, sigmasq=sigma)
#    V<-.Call("kroneckerSpar_shift_OU_EB", R=sig, C=vcvList[[after]], beta=beta, Rrows=as.integer(p),  Crows=as.integer(n), IA=as.integer(IAr), JA=as.integer(JAr), A=as.double(precalcMat@entries))
    
#    loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method="sparse")
    
#       list(loglik=-loglik$logl, theta=loglik$anc)
#}

# sparse matrix eb-ou
lik.shift_ebou_sparse<-function(alpha,sigma,sig,beta,dat,error,vcvList,p){
    Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=vcvList[[before]],lambda=alpha$values,S=alpha$vectors,sigmasq=sigma, S1=alpha$invectors)
    V<-.Call("kroneckerSpar_shift_OU_EB", R=sig, C=vcvList[[after]], beta=matrix(beta,p,p), Rrows=as.integer(p),  Crows=as.integer(n), V=Vou, IA=as.integer(IAr), JA=as.integer(JAr), A=as.double(precalcMat@entries))

    loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc,method="sparse", ch=ch, precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)

       list(loglik=-loglik$logl, theta=loglik$anc)
}

# sparse matrix eb-bm
lik.shift_ebbm_sparse<-function(beta,sigma,sig,dat,error,vcvList,p){
    V<-.Call("kroneckerSpar_shift_EB_BM", R1=sig, R2=sigma, C1=vcvList[[before]], C2=vcvList[[after]], beta=matrix(beta,p,p), Rrows=as.integer(p),  Crows=as.integer(n), IA=as.integer(IAr), JA=as.integer(JAr), A=as.double(precalcMat@entries))
    
    loglik<-loglik_mvmorph(dat,V,W,n,p,error,precalc, method="sparse",ch=ch, precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA)
    
       list(loglik=-loglik$logl, theta=loglik$anc)
}



##----------------------Optimization------------------------------------------##

# select the loglik functions & starting values

switch(model,

"RR"={
    # Brownian and OU models with different rates
    if(method=="sparse"){
    lik.shift<-lik.shift_oubm_sparse
    }else{
    lik.shift<-lik.shift_oubm
    }
    starting<-c(alpha,sigma,sig)
},
"ER"={
    # Brownian and OU models with the same rates
    if(method=="sparse"){
        lik.shift<-lik.shift_oubm_sparse
    }else{
        lik.shift<-lik.shift_oubm
    }
    starting<-c(alpha,sigma)
},
"CV"={
    # Brownian & ACDC models with the same rates
    if(method=="sparse"){
        lik.shift<-lik.shift_ebbm_sparse
    }else{
        lik.shift<-lik.shift_ebbm
    }
    starting<-c(beta,sigma)
},
"CVG"={
    # Brownian & ACDC models with different rates
    if(method=="sparse"){
        lik.shift<-lik.shift_ebbm_sparse
    }else{
        lik.shift<-lik.shift_ebbm
    }
    starting<-c(beta,sigma,sig)
},
"OV"={
    # OU & ACDC models with the same rates
    if(method=="sparse"){
        lik.shift<-lik.shift_ebou_sparse
    }else{
        lik.shift<-lik.shift_ebou
    }
    starting<-c(alpha,sigma,beta)
},
"OVG"={
    # OU & ACDC models with independent rates
    if(method=="sparse"){
        lik.shift<-lik.shift_ebou_sparse
    }else{
        lik.shift<-lik.shift_ebou
    }
    starting<-c(alpha,sigma,sig,beta)
})

# optimizer
if(model=="OVG" || model =="OV"){
     if(optimization!="subplex"){
        estim <- optim(par=starting,fn = function (par) {lik.shift(alpha=buildA(par[seq_len(nalpha)],p),sigma=sym.par(par[nalpha+seq_len(nsigma)]), sig=sym.par(par[nsig+seq_len(nsigma)]), beta=ratevalue(up,low,par[nbeta+1]), error=error, dat=dat, vcvList=vcvList,p=p)$loglik },gr=NULL,hessian=TRUE,method=optimization,control=control)
    }else{
        estim <- subplex(par=starting,fn = function (par) {lik.shift(alpha=buildA(par[seq_len(nalpha)],p),sigma=sym.par(par[nalpha+seq_len(nsigma)]), sig=sym.par(par[nsig+seq_len(nsigma)]), beta=par[nbeta+1], error=error, dat=dat, vcvList=vcvList,p=p)$loglik},hessian=TRUE,control=control)
    }
  
}else if(model=="CV" || model=="CVG"){
    if(optimization!="subplex"){
        estim <- optim(par=starting,fn = function (par) {lik.shift(beta=ratevalue(up,low,par[nbeta]),sigma=sym.par(par[nbeta+seq_len(nsigma)]), sig=sym.par(par[nsig+seq_len(nsigma)]), error=error, dat=dat, vcvList=vcvList,p=p)$loglik},gr=NULL,hessian=TRUE,method=optimization,control=control)
    }else{
        estim <- subplex(par=starting,fn = function (par) {lik.shift(beta=ratevalue(up,low,par[nbeta]),sigma=sym.par(par[nbeta+seq_len(nsigma)]), sig=sym.par(par[nsig+seq_len(nsigma)]), error=error, dat=dat, vcvList=vcvList,p=p)$loglik},hessian=TRUE,control=control)
    }
}else{
     if(optimization!="subplex"){
        estim <- optim(par=starting,fn = function (par) {lik.shift(alpha=buildA(par[seq_len(nalpha)],p),sigma=sym.par(par[nalpha+seq_len(nsigma)]), sig=sym.par(par[nsig+seq_len(nsigma)]), error=error, dat=dat, vcvList=vcvList,p=p)$loglik},gr=NULL,hessian=TRUE,method=optimization,control=control)
    }else{
        estim <- subplex(par=starting,fn = function (par) {lik.shift(alpha=buildA(par[seq_len(nalpha)],p),sigma=sym.par(par[nalpha+seq_len(nsigma)]), sig=sym.par(par[nsig+seq_len(nsigma)]), error=error, dat=dat, vcvList=vcvList,p=p)$loglik},hessian=TRUE,control=control)
    }
}
##---------------------theta estimation---------------------------------------##

if(model=="OVG" || model=="OV"){
res.theta<-lik.shift(alpha=buildA(estim$par[seq_len(nalpha)],p),sigma=sym.par(estim$par[nalpha+seq_len(nsigma)]), sig=sym.par(estim$par[nsig+seq_len(nsigma)]), beta=ratevalue(up,low,estim$par[nbeta+1]), error=error, dat=dat, vcvList=vcvList,p=p)$theta
}else if(model=="CV" || model=="CVG"){
res.theta<-lik.shift(beta=ratevalue(up,low,estim$par[nbeta]),sigma=sym.par(estim$par[nbeta+seq_len(nsigma)]), sig=sym.par(estim$par[nsig+seq_len(nsigma)]), error=error, dat=dat, vcvList=vcvList,p=p)$theta
}else{
res.theta<-lik.shift(alpha=buildA(estim$par[seq_len(nalpha)],p),sigma=sym.par(estim$par[nalpha+seq_len(nsigma)]), sig=sym.par(estim$par[nsig+seq_len(nsigma)]), error=error, dat=dat, vcvList=vcvList,p=p)$theta
}
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
# Loglikelihood estimate
LL<- -estim$value

# free parameters
nparam=p+length(estim$par)

# maximum likelihood estimates of alpha and sigma
estim.alpha<-estim$par[seq_len(nalpha)]
estim.sigma<-estim$par[nalpha+seq_len(nsigma)]
if(model=="RR" || model=="CVG" || model=="OVG"){
estim.sig<-estim$par[nsig+seq_len(nsigma)]
}


# alpha matrix
alpha.mat<-buildA(estim.alpha,p)$A
rownames(alpha.mat)<-colnames(data)
colnames(alpha.mat)<-colnames(data)
# sigma matrix
sigma.mat<-sym.par(estim.sigma)
rownames(sigma.mat)<-colnames(data)
colnames(sigma.mat)<-colnames(data)
# sig matrix
if(model=="RR" || model=="CVG" || model=="OVG"){
sig.mat<-sym.par(estim.sig)
rownames(sig.mat)<-colnames(data)
colnames(sig.mat)<-colnames(data)
}
# beta matrix
if(model=="OVG" || model=="OV"){
estim.beta<-ratevalue(up,low,estim$par[nbeta+1])
beta.mat<-c(estim.beta)
}
if(model=="CV" || model=="CVG"){
estim.beta<-ratevalue(up,low,estim$par[nbeta])
beta.mat<-c(estim.beta)
}

# AIC
AIC<- -2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1989
# matrix of estimated theta values
theta.mat<-matrix(res.theta,1)

rownames(theta.mat)<-"theta"
colnames(theta.mat)<-colnames(data)

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("-- Summary results for the",nmod," --","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters")
cat("\n")
cat("Estimated theta values","\n")
cat("______________________","\n")
print(theta.mat)
cat("\n")
if(model=="CV" || model=="CVG"){
    cat("ML beta values","\n")
    cat("______________________","\n")
    print(beta.mat)
}else{
    cat("ML alpha values","\n")
    cat("______________________","\n")
    print(alpha.mat)
}
cat("\n")
cat("ML sigma values","\n")
cat("______________________","\n")
print(sigma.mat)
if(model=="RR" || model=="radiate"){
    cat("\n")
    cat("ML sigma radiation values","\n")
    cat("______________________","\n")
    print(sig.mat)
}
if(model=="CVG" || model=="OVG"){
    cat("\n")
    cat("ML sigma values ( recent slice:",colnames(tree$mapped.edge)[2],")","\n")
    cat("______________________","\n")
    print(sig.mat)
}
if(model=="OVG" || model=="OV"){
    cat("\n")
    cat("ML beta values","\n")
    cat("______________________","\n")
    print(beta.mat)
}
}

##-------------------Save infos in parameters---------------------------------##
param$nparam<-nparam
param$nbspecies<-n
param$ntraits<-p
param$nregimes<-2
param$method<-method
param$optimization<-optimization
param$traits<-colnames(data)
param$names_regimes<-colnames(tree$mapped.edge)
param$model<-c(model,nmod)
param$before<-before
param$after<-after
##------------------List results----------------------------------------------##
if(model=="ER"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha=alpha.mat, sigma=sigma.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
}else if(model=="RR"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha=alpha.mat, sigma=sigma.mat, sig=sig.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
}else if(model=="CV"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, beta=alpha.mat, sigma=sigma.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
}else if(model=="CVG"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, beta=alpha.mat, sigma=sigma.mat, sig=sig.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
}else if(model=="OV"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha=alpha.mat, sigma=sigma.mat, beta=beta.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
}else if(model=="OVG"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha=alpha.mat, sigma=sigma.mat, sig=sig.mat, beta=beta.mat, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val, param=param)
}
class(results)<-c("mvmorph","shift")
invisible(results)

    }
} #fin de la fonction
