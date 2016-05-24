################################################################################
##                                                                            ##
##                               mvMORPH: fun.r                               ##
##                                                                            ##
##   Internal functions for the mvMORPH package                               ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013 - updated 10-03-2015                ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam                           ##
##                                                                            ##
################################################################################

##------------------------Fonctions necessaires-------------------------------##

# Calcul d'une matrice positive semi-definie (Choleski transform) modified from OUCH by A. King
sym.par<-function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    stop("a symmetric matrix is parameterized by a triangular number of parameters",call.=FALSE)
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  tcrossprod(y)
}

# Inverse of Sym.par, return a vector (from OUCH) by A. King
sym.unpar <- function (x) {
  y <- t(chol(x))
  y[lower.tri(y,diag=TRUE)]
}

sym.unpar_off <- function (x) {
    y <- t(chol(x))
    y[lower.tri(y,diag=FALSE)]
}



# alpha matrix parameterization
matrixParam<-function(x,p,matrix="symmetric",tol=0.000001){
    
    switch(matrix,
    "symmetricPositive"={ # Cholesky decomposition
        y <- matrix(0,p,p)
        y[lower.tri(y,diag=TRUE)] <- x
        A<-tcrossprod(y)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "symmetric"={ # SVD decomposition
        dim1<-p*(p-1)/2
        Q<-.Call("givens_ortho", Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-x[(dim1+1):(dim1+p)]+tol
        invQ<-t(Q)
        A<-Q%*%diag(T)%*%invQ
        Adecomp<-list(vectors=Q, values=T, A=A, invectors=invQ)
    },
    "nsymmetric"={ # Schur decomposition
        dim1<-p*(p-1)/2
        Q<-.Call("givens_ortho", Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-diag(x[(dim1+1):(dim1+p)],p)
        T[upper.tri(T)]<-x[(dim1+p+1):(p*p)]
        A<-Q%*%T%*%t(Q)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "nsymPositive"={ # Schur decomposition
        dim1<-p*(p-1)/2
        Q<-.Call("givens_ortho", Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-diag(exp(x[(dim1+1):(dim1+p)]),p) # exp to force positive eigenvalues
        T[upper.tri(T)]<-x[(dim1+p+1):(p*p)]
        A<-Q%*%T%*%t(Q)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "diagonal"={
        A<-diag(x,p)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "diagonalPositive"={
        A<-diag(exp(x)+tol,p)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "lower"={
        dim1<-p*(p-1)/2
        A<-matrix(0,p,p)
        A[lower.tri(A,diag=FALSE)] <- x[1:dim1]
        diag(A) <- exp(x[(dim1+1):(dim1+p)])
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "upper"={
        dim1<-p*(p-1)/2
        A<-matrix(0,p,p)
        A[upper.tri(A,diag=FALSE)] <- x[1:dim1]
        diag(A) <- exp(x[(dim1+1):(dim1+p)])
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "univariate"={
        Adecomp<-list(vectors=1, values=x, A=x, invectors=1)
    }
    )
    return(Adecomp)
}


# Compute factor list / Creation d'une liste de facteurs
newList<-function(factorVal,nv){
    temp<-factorVal
    factorVal<-lapply(1:(nv-1),function(i){val<-c(max(factorVal)+temp);val<-c(factorVal,val);val})
    factorVal<-unlist(factorVal)
    return(factorVal)
}

# Compute the design matrix for multiple mean / Creation de la matrice de variables indicatrices pour plusieurs moyennes

multD<-function(tree,k,nbtip,smean=TRUE){
    if(smean==TRUE){
        y<-matrix(0, nbtip * k, k)
        for (i in 1:(nbtip * k)){
            for (j in 1:k){
                if ((j - 1) * nbtip < i && i <= j * nbtip){y[i, j] = 1 }
            }
        }
    }else{
        if(is.null(tree[["mapped.edge"]])){
          stop("The specified tree must be in SIMMAP format with different regimes")
        }
        namestip<-sapply(1:nbtip,function(i){
        ind<-which(tree$edge[,2]==i);
        names(tree$maps[[ind]][length(tree$maps[[ind]])])})
        group<-as.numeric(as.factor(namestip))
        
        if(k==1){
            ngr<-length(unique(group))
            y=matrix(0,ncol=ngr,nrow=nbtip)
            for (i in 1:nbtip){y[i,group[i]] <- 1}
        }else{
        gr<-newList(group,k)
        ngr<-unique(gr)
        tgr<-length(gr)
        y=matrix(0,tgr,length(ngr))
        for (i in 1:tgr){y[i,gr[i]] <- 1}
        }
    }
return(y)
}


# Compute matrix time to common ancestor
mTime<-function(phy,scale.height){
  if(is.ultrametric(phy)){
   if(scale.height==TRUE){
   mdist<-vcv.phylo(phy)/max(vcv.phylo(phy))
   }else{
    mdist<-vcv.phylo(phy)}
    mSpdist<-rep(max(mdist),length(phy$tip.label))
  }else{
   vcv<-vcv.phylo(phy)
   if(scale.height==TRUE){
   vstand<-vcv/max(vcv)
   }else{vstand<-vcv}
   mSpdist<-diag(vstand)
   mdist<-vstand
  }
  # mcoph<-diag(mdist)-mdist
list(mSpDist=mSpdist, mDist=mdist)
}

# Binary regime coding
regimeList<-function(mm,k,root=TRUE){ # use root="stationary" for ouch...
    nReg=length(mm)
    regime <- matrix(0,nrow=nReg,ncol=k) # remplacer 0?
    for(i in 1:nReg){
        regime[i,mm[i]]<-1
    }
    if(root=="stationary"){
        regime[nReg,]<-0 ## under OU1 implicitely assumed stationary, just use root=TRUE or FALSE (ultrametric trees)
    }
    return(regime)
}

# Set root regime
indiceReg<-function(n,indice,facInd, root=TRUE){
    if(root==TRUE){
        for(i in 1:n){
            val=length(indice[[i]])
            indice[[i]][val+1]<-facInd[facInd=="_root_state"]
        }
    }else{
        for(i in 1:n){
            val=length(indice[[i]])
            indice[[i]][val+1]<-indice[[i]][val]
        }
    }
    return(indice)
}


# Test for polytomies
eval_polytom<-function(tree){
    nb.tip <- length(tree$tip.label)
    nb.node <- tree$Nnode
    if (nb.node != nb.tip - 1) {
        stop("You can't use \"pic\" with polytomies, try instead \"rpf\",\"inverse\",\"pseudoinverse\" or \"sparse\". Otherwise, first transform the tree using the \"multi2di\" function")
    }
}

# Estimate starting values for the variance matrix
varBM<-function(tree,data,n,k){
    if(any(is.na(data))){
        res<-diag(0.1,k)
        return(res)
    }
    nb.tip <- length(tree$tip.label)
    nb.node <- tree$Nnode
    if (nb.node != nb.tip - 1) {
    tree <- multi2di(tree)
    }
    rate<-rep(0,k)
    tree<-reorder.phylo(tree,"postorder")
    value<-list(tree$edge.length)
    res<-.Call("PIC_gen", x=as.vector(as.matrix(data)), n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=value, times=1, rate=rate, Tmax=1, Model=as.integer(10), mu=1, sigma=1)
    return(res[[2]])
}

# Generate random multivariate distributions
rmvnorm<-function(n=1, mean, var){
    p<-length(mean)
    if (!all(dim(var)==c(p,p)))
    stop("length of ",sQuote("mean")," must equal the dimension of the square matrix ",sQuote("var"))
    cf <- t(chol(var))
    matrix(mean,p,n)+cf%*%matrix(rnorm(p*n),p,n)
}

# Generate a multi-phylo list for SIMMAP trees
vcvSplit<-function(tree){
    multi.tre<-list()
    class(multi.tre)<-"multiPhylo"
    #Array method - better for memory?
    #C2<-array(dim=c(nrow(C1),ncol(C1),ncol(tree$mapped.edge)))
    C<-list()
    for(i in 1:ncol(tree$mapped.edge)){
        multi.tre[[i]]<-tree
        multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
        multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
        temp<-vcv.phylo(multi.tre[[i]])
        C[[i]]<-temp
    }
return(C)
}

# Generate box constraints for the likelihood search
ratevalue<-function(up,low,x){
    y<-x<low | x>up
    x[y]=0
    return(x)
}

# Compute the stationary multivariate normal distribution for the multivariate Ornstein-Uhlenbeck (Bartoszek et al. 2012 - B.8)
StationaryVariance <- function(alpha,sigma){
    sigma <- sigma
    eig <- eigen(alpha)
    P <- eig$vectors
    invP <- solve(P)
    eigvalues <- eig$values
    p=dim(sigma)[1]
    Mat <- matrix(0,p,p)
    for(i in 1:p){
        for(j in 1:p){
            Mat[i,j] <- 1/(eigvalues[i]+eigvalues[j])
        }
    }
    StVar <- P%*%(Mat*(invP%*%sigma%*%t(invP)))%*%t(P)
    return(StVar)
}

##-----------------Function used in mvBM--------------------------------------##


# Constrainded Choleski decomposition (Adams, 2013: Systematic Biology)                   
build.chol<-function(b,p){
 c.mat<-matrix(0,nrow=p,ncol=p)
 c.mat[lower.tri(c.mat)] <- b[-1]
 c.mat[p,p]<-exp(b[1])
 c.mat[1,1]<-sqrt(sum((c.mat[p,])^2))
 if(p>2){
 for (i in 2:(p-1)){
 c.mat[i,i]<-ifelse((c.mat[1,1]^2-sum((c.mat[i,])^2))>0,sqrt(c.mat[1,1]^2-sum((c.mat[i,])^2)), 0)
 }}
 return(c.mat)
}


##----------------------mvfit_likelihood--------------------------------------##

loglik_mvmorph<-function(data,V=NULL,D=NULL,n,k,error=NULL,precalc=precalc,method, param=list(),ch=NULL,precalcMat=NULL,sizeD=NULL,NA_val=NULL,Indice_NA=NULL){

data<-as.numeric(as.matrix(data))
size<-k*n
if(NA_val==TRUE){
    V<-V[-Indice_NA,-Indice_NA]
    D<-D[-Indice_NA,]
    data<-data[-Indice_NA]
    size<-length(data)
}

switch(method,

"rpf"={
    if(is.null(error)!=TRUE){ ms<-1 }else{ ms<-0} # length of the error vector should be the same as the data vector
    cholres<-.Call("Chol_RPF",V,D,data,as.integer(sizeD),as.integer(size),mserr=error,ismserr=as.integer(ms))
    beta<-pseudoinverse(cholres[[3]])%*%cholres[[4]]
    det<-cholres[[2]]
    residus=D%*%beta-data
    quad<-.Call("Chol_RPF_quadprod", cholres[[1]], residus, as.integer(size))
    logl<--.5*quad-.5*as.numeric(det)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=beta)
},
"univarpf"={
    if(is.null(error)!=TRUE){ ms<-1 }else{ ms<-0}
    size<-k*n
    cholres<-.Call("Chol_RPF_univ",V,D,data,as.integer(sizeD),as.integer(size),mserr=error,ismserr=as.integer(ms))
    beta<-pseudoinverse(cholres[[3]])%*%cholres[[4]]
    det<-cholres[[2]]
    residus=D%*%beta-data
    quad<-.Call("Chol_RPF_quadprod_column", cholres[[1]], residus, as.integer(size))
    logl<--.5*quad-.5*as.numeric(det)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=beta)
},
"sparse"={
    ## On considere que les valeurs dans @entries de V de precalc ont ete modifiees / assume that the values in @entries were updated
    if(is.null(error)==FALSE){
        if(!is.null(precalc)){
        diag(precalc$V)<-diag(precalc$V)+error
        }else{
        diag(precalcMat)<-diag(precalcMat)+error
        }
    }
    if(!is.null(precalc)){
    U<-update.spam.chol.NgPeyton(precalc$ch,precalc$V)
    }else{
    U<-update.spam.chol.NgPeyton(ch,precalcMat)
    }
    vec<-forwardsolve(U,data)
    xx<-forwardsolve(U,D)
    beta<-pseudoinverse(matrix(xx,ncol=sizeD))%*%vec
    res<-D%*%beta-data
    vec1<-forwardsolve(U,res)
    a<-sum(vec1^2)
    DET<-determinant(U)
    logl<--.5*a-.5*as.numeric(DET$modulus*2)-.5*(n*k*log(2*pi))
    results<-list(logl=logl,anc=beta)
},
"pseudoinverse"={
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
    
    inv<-pseudoinverse(V)
    beta<-pseudoinverse(t(D)%*%inv%*%D)%*%t(D)%*%inv%*%data
    DET<-determinant(V, logarithm=TRUE)
    res<-D%*%beta-data
    logl<--.5*(t(res)%*%inv%*%(res))-.5*as.numeric(DET$modulus)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=beta)
},
"inverse"={
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }

    inv<-solve(V)
    beta<-solve(t(D)%*%inv%*%D)%*%t(D)%*%inv%*%data
    DET<-determinant(V, logarithm=TRUE)
    res<-D%*%beta-data
    logl<--.5*(t(res)%*%inv%*%(res))-.5*as.numeric(DET$modulus)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=beta)
})


return(results)


}


##----------------------Print_functions---------------------------------------##

print.bm<-function(x,...){
    
    cat("\n")
    if(x$param$constraint==TRUE){
        cat("-- Summary results for multiple constrained rates",x$param$model,"model --","\n")
    }else if(x$param$constraint=="correlation"){
        cat("-- Summary results for common correlation ",x$param$model,"model --","\n")
    }else if(x$param$constraint=="shared"){
        cat("-- Summary results for shared eigenvectors ",x$param$model,"model --","\n")
    }else if(x$param$constraint=="proportional"){
        cat("-- Summary results for proportional rates matrices ",x$param$model,"model --","\n")
    }else{
        cat("-- Summary results for multiple rates",x$param$model,"model --","\n")
    }
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters","\n")
    cat("\n")
    cat("Estimated rates matrix","\n")
    cat("______________________","\n")
    print(x$sigma)
    cat("\n")
    cat("Estimated root state","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")

}

print.acdc<-function(x,...){
    cat("\n")
    cat("-- Summary results for Early Burst or ACDC model --","\n")
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters","\n")
    cat("Rate change:","\n")
    cat("______________________","\n")
    print(x$beta)
    cat("\n")
    cat("Estimated rates matrix","\n")
    cat("______________________","\n")
    print(x$sigma)
    cat("\n")
    cat("Estimated root states","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
}

print.ou<-function(x,...){
    cat("\n")
    cat("-- Summary results --","\n")
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters","\n")
    cat("\n")
    cat("Estimated theta values","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
    cat("ML alpha values","\n")
    cat("______________________","\n")
    print(x$alpha)
    cat("\n")
    cat("ML sigma values","\n")
    cat("______________________","\n")
    print(x$sigma)
}

print.shift<-function(x,...){
    cat("\n")
    cat("-- Summary results for the",x$param$model[2]," --","\n")
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters")
    cat("\n")
    cat("Estimated theta values","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
    if(x$param$model[1]=="CV" || x$param$model[1]=="CVG"){
        cat("ML beta values","\n")
        cat("______________________","\n")
        print(x$beta)
    }else{
        cat("ML alpha values","\n")
        cat("______________________","\n")
        print(x$alpha)
    }
    cat("\n")
    cat("ML sigma values","\n")
    cat("______________________","\n")
    print(x$sigma)
    if(x$param$model[1]=="RR" || x$param$model[1]=="radiate"){
        cat("\n")
        cat("ML sigma radiation values","\n")
        cat("______________________","\n")
        print(x$sig)
    }
    if(x$param$model[1]=="CVG" || x$param$model[1]=="OVG"){
        cat("\n")
        cat("ML sigma values ( recent slice:",x$param$names_regimes[2],")","\n")
        cat("______________________","\n")
        print(x$sigma)
    }
    if(x$param$model[1]=="OVG" || x$param$model[1]=="OV"){
        cat("\n")
        cat("ML beta values","\n")
        cat("______________________","\n")
        print(x$beta)
    }
}

print.mvmorph.lrt<-function(x,...){
    if(x$pval<0.000001){signif<-c("***")}else if(x$pval<0.001){
        signif<-c("**") }else if(x$pval<0.01){signif<-c("*")}else if(x$pval<0.05){signif<-c(".")}else{signif<-""}
    cat("-- Log-likelihood Ratio Test --","\n")
    cat("Model",x$model1," versus ",x$model2,"\n")
    cat("Number of degrees of freedom :",x$ddf,"\n")
    cat("LRT statistic:",x$ratio," p-value:",x$pval,signif,"\n")
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}

print.mvmorph.precalc<-function(x,...){
    cat("A tree with",length(x$tree$tip.label),"species used in precalc","\n")
    cat("Optimized for:","\n")
    cat("-----------------","\n")
    cat("method:",x$param$method,"\n")
    cat("model:",x$model,"\n")
    cat("number of traits:",x$param$nbtraits,"\n")
}

summary.mvmorph<-function(object,...){
    cat("mvMORPH model :",object$param$model," summary","\n")
    cat("AIC :",object$AIC,"\n")
    cat("AICc:",object$AICc,"\n")
    cat("Log-Likelihood:",object$LogLik,"\n")
    if(object$convergence==0){cat("Succesful convergence","\n")}else{cat("Convergence has not been reached","\n")}
    if(object$hess.value==0){cat("Reliable solution","\n")}else{cat("Unreliable solution (Likelihood at a saddle point)","\n")}
}

## Return the model AIC
AIC.mvmorph<-function(object,...,k){
    return(object$AIC)
}

## Return the model AICc
AICc <- function(x) UseMethod("AICc")
AICc.mvmorph<-function(object){
    return(object$AICc)
}

## Return the model fit log-likelihood
logLik.mvmorph<-function(object,...){
    return(object$LogLik)
}

## Change to include the tree in the analysis? == problematic for large trees... need too much storage!
simulate.mvmorph<-function(object,nsim=1,seed=NULL,...){
    mvSIM(...,param=object,nsim=nsim)
}

## Return the stationary variance for the multivariate Ornstein-Uhlenbeck

stationary.mvmorph<-function(object){
    if(any(class(object)=="shift") | any(class(object)=="ou") ){
            if(is.null(object[["alpha"]])==TRUE){
                stop("The stationary distribution can be computed only for models including Ornstein-Uhlenbeck processes.")
            }
    statMat<-StationaryVariance(object$alpha,object$sigma)
    rownames(statMat)<-rownames(object$sigma)
    colnames(statMat)<-colnames(object$sigma)
    return(statMat)   
    }else{
      warning("The stationary distribution can be computed only for Ornstein-Uhlenbeck processes.","\n")
    }
}

## Compute the phylogenetic half-life

halflife.mvmorph<-function(object){
    if(class(object)[2]=="ou"){
    lambda<-eigen(object$alpha)$values
    phyhalflife<-log(2)/lambda
        return(phyhalflife)
    }else{
        warning("The phylogenetic half-life is computed only for Ornstein-Uhlenbeck models.","\n")  
    }
}