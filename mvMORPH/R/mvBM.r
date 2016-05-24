################################################################################
##                                                                            ##
##                               mvMORPH: mvBM                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################


mvBM<-function(tree, data, error=NULL, model=c("BMM","BM1"),param=list(constraint=FALSE, smean=TRUE), method=c("rpf","pic","sparse","inverse","pseudoinverse"), scale.height=FALSE, optimization=c("L-BFGS-B","Nelder-Mead","subplex"), control=list(maxit=20000), precalc=NULL, diagnostic=TRUE, echo=TRUE){

# select default model
model<-model[1]
method<-method[1]
#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
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

if(is.null(param[["constraint"]])==TRUE){constraint<-param$constraint<-FALSE}else{
    if((param$constraint=="shared" & model=="BM1") | (param$constraint=="proportional" & model=="BM1") | (param$constraint=="correlation" & model=="BM1")){
    constraint<-FALSE
    if(echo==TRUE) cat(" \"shared\" and \"proportional\" can be used only with BMM model","\n")

    }else if(param$constraint=="shared" & k==1){
        constraint<-FALSE
    }else{
    constraint<-param$constraint
    }
}


##------------------------Precalc---------------------------------------------##
if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
    tree<-precalc$tree
    
    if(is.null(tree[["mapped.edge"]])==TRUE){
        model<-"BM1"
       if(echo==TRUE) cat("No selective regimes mapped on the tree, only a BM1 model could be estimated","\n")
    }
    C1<-precalc$C1
    D<-precalc$D
    if(model=="BM1"){param$smean<-TRUE}
    if(model=="BMM"){
    C2<-precalc$C2
    }
    # number of selective regimes
    if(model!="BM1"){
        p<-length(C2)
    }else{ p<-1 }
    
    if(method=="sparse"){
        # Yale sparse format
        JAr<-precalc$JAr
        IAr<-precalc$IAr
        ch<-precalc$ch
        precalcMat<-precalc$V
    }

    
}else{
##------------------------Precalc-off-----------------------------------------##
##------------------------Create VCV matrix-----------------------------------##
if(is.null(tree[["mapped.edge"]])==TRUE){
    model<-"BM1"
    if(echo==TRUE) cat("No selective regimes mapped on the tree, only a BM1 model could be estimated","\n")
}
# Scale the tree
if(scale.height==TRUE){
    maxHeight<-max(nodeHeights(tree))
    tree$edge.length<-tree$edge.length/maxHeight
    if(model=="BMM")tree$mapped.edge<-tree$mapped.edge/maxHeight
}

if(method=="pic"){
    ind<-reorder.phylo(tree,"postorder", index.only=TRUE)
    tree$edge<-tree$edge[ind,]
    tree$edge.length<-tree$edge.length[ind]
    value<-list(tree$edge.length)
    
    if(model=="BMM"){
        tree$mapped.edge<-tree$mapped.edge[ind,]
    }
    # method for computing the log-likelihood
    if(model=="BMM" & k!=1){
        #mvMORPH-1.0.3
        warning("Sorry, the \"pic\" method only works with univariate data for the BMM model ","\n","the \"rpf\" method has been used instead...","\n")
        method<-"rpf"
    }
    C1<-NULL
    C2<-NULL
}

if(method!="pic"){
# Compute vcv for SIMMAP tree (replace in precalc?)
    C1<-vcv.phylo(tree)
	if(!is.null(rownames(data))) {
	 if(any(tree$tip.label==rownames(data))){
         C1<-C1[rownames(data),rownames(data)]
  }else if(echo==TRUE){
    cat("row names of the data matrix must match tip names of your phylogeny!","\n")
    }}else if(echo==TRUE){
	cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
	}
if(model=="BMM"){
  multi.tre<-list()
  class(multi.tre)<-"multiPhylo"
  #Array method
  #C2<-array(dim=c(nrow(C1),ncol(C1),ncol(tree$mapped.edge)))
  C2<-list()
	for(i in 1:ncol(tree$mapped.edge)){
		multi.tre[[i]]<-tree
		multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		if(any(tree$tip.label==rownames(data))) { 
		C2[[i]]<-temp[rownames(data),rownames(data)]
		}else{
            C2[[i]]<-temp
        }
      }
    }
        
}

##------------------------Parameters------------------------------------------##

# number of selective regimes
if(model!="BM1"){
    p<-length(colnames(tree$mapped.edge))
}else{ p<-1 }

        
# Compute the design matrix
if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
        if(model=="BM1"){param$smean<-TRUE}
        D<-multD(tree,k,n,smean=param$smean)

if(method=="sparse"){
    V<-kronecker((matrix(1,k,k)+diag(k)), C1)
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


}#end off-precalc

# Define the variance-covariance functions
if(model=="BMM"){
switch(method,
"rpf"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kroneckerSum", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n), dimlist=as.integer(p))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA)
    return(loglik)
    }
},
"sparse"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kroneckerSumSpar", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n),  dimlist=as.integer(p), IA=IAr, JA=JAr, A=precalcMat@entries)
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=FALSE,Indice_NA=NULL)
    return(loglik)
    }
},
"pseudoinverse"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kroneckerSum", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n), dimlist=as.integer(p))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA)
    return(loglik)
    }
},
"inverse"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kroneckerSum", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n), dimlist=as.integer(p))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA)
    return(loglik)
    }
},
"pic"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    # Mult only available for univariate case
    k=1
    sig<-unlist(sig)
    tree$edge.length<-tree$mapped.edge%*%sig
    # Compute the LLik
    res<-.Call("PIC_gen", x=dat, n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=list(tree$edge.length), times=1, rate=rep(0,k), Tmax=1, Model=as.integer(7), mu=1, sigma=1)
    logl<- -0.5 * ( n * k * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
    return(list(logl=logl,ancstate=res[[7]], sigma=res[[2]]))
        }
})

}else{

switch(method,
"rpf"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kronecker_mvmorph", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA)
    return(loglik)
    }
},
"sparse"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kroneckerSumSpar", R=list(sig), C=list(C), Rrows=as.integer(k),  Crows=as.integer(n),  dimlist=as.integer(1), IA=IAr, JA=JAr, A=precalcMat@entries)
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=FALSE,Indice_NA=NULL)
    return(loglik)
    }
},
"pseudoinverse"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kronecker_mvmorph", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA)
    return(loglik)
    }
},
"inverse"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    V<-.Call("kronecker_mvmorph", R=sig, C=C, Rrows=as.integer(k),  Crows=as.integer(n))
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA)
    return(loglik)
    }
},
"pic"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,k,error,method){
    res<-.Call("PIC_gen", x=dat, n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=list(tree$edge.length), times=1, rate=rep(0,k), Tmax=1, Model=as.integer(7), mu=NULL, sigma=sig)
    logl<- -0.5 * ( n * k * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
    return(list(logl=logl,ancstate=res[[7]], sigma=res[[2]]))
    }
}
)
# End if
}

# sigma matrix prameterization

buildSigma<-function(par,index.mat=NULL,sig=NULL,model,constraint){
    if(constraint==FALSE){
        switch(model,
        "BMM"={
            sig[]<-c(par)[index.mat]
            sigma<-lapply(1:p,function(x){ sym.par(sig[x,])})
        },
        "BM1"={
            sigma<-sym.par(par)
        })
    }else if(constraint=="diagonal"){
        
        switch(model,
        "BMM"={
            sig[] <- c(par)[index.mat]
            sigma<-lapply(1:p,function(x){ diag(diag(sig[x,]%*%t(sig[x,])))})
        },
        "BM1"={
            sigma<-diag(diag(par%*%t(par)))
        })
        
    }else if(constraint=="shared"){
        # Comparisons between groups (ie not allowed for BM1)
        
        diagval<-1:k
        variance <- par[diagval]
        param<-par[-diagval]
        sig[] <- c(param)[index.mat]
        sigma<-lapply(1:p,function(x){.Call("spherical", param=sig[x,], variance=variance, dim=as.integer(k))})
        
    }else if(constraint=="correlation"){
        # Comparisons between groups (ie not allowed for BM1)
        
        diagval<-1:npar
        angle <- par[diagval]
        param<-par[-diagval]
        sig[] <- c(param)[index.mat]
        sigma<-lapply(1:p,function(x){.Call("spherical", param=angle, variance=sig[x,], dim=as.integer(k))})
        
    }else if(constraint=="proportional"){
        # Comparisons between groups (ie not allowed for BM1)
        propindex<-1:(p-1)
        constant<-par[propindex]
        rates_mat<-sym.par(par[-propindex])
        sigma<-lapply(1:p,function(x){ if(x==1){rates_mat}else{rates_mat*as.vector((constant[x-1]%*%t(constant[x-1])))} })
      
    }else{
        switch(model,
        "BMM"={
            sig[] <- c(par)[index.mat]
            sigma<-lapply(1:p,function(x){ tcrossprod(build.chol(sig[x,],k))})
        },
        "BM1"={
            sigma<-tcrossprod(build.chol(par,k))
        })
    }
    
    return(sigma)
}


##------------------LogLikelihood function for multiple rates per traits------##
lik.Mult<-function(par,dat,C,D,index.mat,sig,error, p, k, n, precalcMat, method,constraint){

  loglik<-bm_fun_matrix(C,buildSigma(par,index.mat,sig,model,constraint),dat,D,precalcMat,n,k,error,method)
  
  list(loglik=-loglik$logl, ancstate=loglik$anc)
  
}
##---------------------Loglik BM1---------------------------------------------##
lik.BM1<-function(par,dat,C,D,error,method,precalcMat,n,k, constraint){ ##
  
  loglik<-bm_fun_matrix(C,buildSigma(par,NULL,NULL,model,constraint),dat,D,precalcMat,n,k,error,method)
  
  list(loglik=-loglik$logl, ancstate=loglik$anc)
}


if(model=="BMM"){
        if(param$constraint==FALSE){
##---------------------Optimization BMM---------------------------------------##
        # number of parameters
        npar=(k*(k+1)/2)
        # sigma matrix
        sig<-matrix(1,p,npar)

        # index matrix of rates
        index.mat<-matrix(1:length(sig),p,npar,byrow=TRUE)

        # initial values for the optimizer
        if(is.null(param[["sigma"]])==TRUE){
            sig1<-varBM(tree,data,n,k)
            sig1<-sym.unpar(sig1)
            starting<-unlist(lapply(1:p,function(x){sig1}))
        }else{
            if(length(param$sigma[[1]])==npar){
                starting<-unlist(param$sigma)
            }else{
                starting<-unlist(lapply(1:length(param$sigma),function(x){sym.unpar(param$sigma[[x]])}))
            }
            if(length(starting)!=(p*npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
        }
        
##---------------------Optimization BMM diagonal------------------------------##
        }else if(param$constraint=="diagonal"){
            # number of parameters
            npar=k
            # sigma matrix
            sig<-matrix(1,p,npar)
            
            # index matrix of rates
            index.mat<-matrix(1:length(sig),p,npar,byrow=TRUE)
            
            # initial values for the optimizer
            if(is.null(param[["sigma"]])==TRUE){
                sig1<-varBM(tree,data,n,k)
                sig1<-diag(sig1)
                starting<-unlist(lapply(1:p,function(x){sig1}))
            }else{
                if(length(param$sigma[[1]])==npar){
                    starting<-unlist(param$sigma)
                }else{
                    starting<-unlist(lapply(1:length(param$sigma),function(x){diag(param$sigma[[x]])}))
                }
                if(length(starting)!=(p*npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
            }
##---------------------Optimization BMM shared eigenvectors-------------------##
        }else if(param$constraint=="shared"){
            
            # number of parameters
            npar=(k*(k-1)/2)
            # sigma matrix
            sig<-matrix(1,p,npar)
            
            # index matrix of rates
            index.mat<-matrix(1:length(sig),p,npar,byrow=TRUE)

            # initial values for the optimizer
            if(is.null(param[["sigma"]])==TRUE){
                sig1<-varBM(tree,data,n,k)
                varval<-diag(sig1)
                sig1<-sym.unpar_off(sig1)
                starting<-unlist(lapply(1:p,function(x){sig1}))
                starting<-c(varval,starting)
            
            }else{
                if(length(param$sigma[[1]])==npar+k){
                    varval<-diag(sym.par(param$sigma[[1]]))
                    starting<-lapply(1:length(param$sigma),function(x){sym.unpar_off(sym.par(param$sigma[[x]]))})
                }else{
                   starting<-unlist(lapply(1:length(param$sigma),function(x){sym.unpar_off(param$sigma[[x]])}))
                   varval<-diag(param$sigma[[1]])
                }
               starting<-c(varval,starting)
                if(length(starting)!=p*npar+k){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
            }
##---------------------Optimization BMM shared eigenvectors-------------------##
        }else if(param$constraint=="correlation"){
            
            # number of parameters
                npar=(k*(k-1)/2)
            # sigma matrix
                sig<-matrix(1,p,k)
            
            # index matrix of rates
                index.mat<-matrix(1:length(sig),p,k,byrow=TRUE)
            
            # initial values for the optimizer
             if(is.null(param[["sigma"]])==TRUE){
                 sig1<-varBM(tree,data,n,k)
                 starting<-sym.unpar_off(sig1)
                 varval<-unlist(lapply(1:p,function(x){diag(sig1)}))
                 starting<-c(starting,varval)
             
                }else{
                if(length(param$sigma[[1]])==npar+k){
                    starting<-sym.unpar_off(sym.par(param$sigma[[1]]))
                    varval<-lapply(1:length(param$sigma),function(x){diag(sym.par(param$sigma[[1]]))})
                }else{
                    starting<-sym.unpar_off(param$sigma[[1]])
                    varval<-unlist(lapply(1:length(param$sigma),function(x){diag(param$sigma[[x]])}))
                }
                starting<-c(starting,varval)
                if(length(starting)!=p*k+npar){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
                }
            

            
##---------------------Optimization BMM proportional--------------------------##
        }else if(param$constraint=="proportional"){
            
              npar=(k*(k+1)/2)
            # initial values for the optimizer
             if(is.null(param[["sigma"]])==TRUE){
                sig1<-varBM(tree,data,n,k)
                starting<-sym.unpar(sig1)
            }else{
                if(length(param$sigma)==npar){
                    starting<-param$sigma
                }else{
                    starting<-sym.unpar(param$sigma)
                }
                if(length(starting)!=(npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
            }
            proportionalval<-rep(1,p-1)
            starting<-c(proportionalval,starting)
            index.mat<-NULL
       
        }else{
##---------------------Optimization BMM constrained---------------------------##
        # number of parameters for the constrained model
        npar=(k*(k-1)/2)+1
        # sigma matrix
        sig<-matrix(1,p,npar)
        # index matrix
        index.mat<-matrix(1:length(sig),p,npar,byrow=TRUE)
        
        if(is.null(param[["sigma"]])==TRUE){
            # initial values for the optimizer
            sig1<-varBM(tree,data,n,k)
            # starting values following Adams (2012)
            sigma.mn<-mean(diag(sig1))
            R.offd<-rep(0,(k*(k-1)/2))
            # # m?me chose mais pour chaque regimes dans la matrice index
            valstart=c(sigma.mn,R.offd)
            starting<-NULL
            for(i in 1:p){
                starting<-c(starting,valstart)
            }
        }else{ ## starting values are provided
            if(length(param$sigma[[1]])==npar){
                starting<-unlist(param$sigma)
            }else{
                starting<-unlist(lapply(1:length(param$sigma),function(x){c(param$sigma[[x]][[1]],param$sigma[[x]][lower.tri(param$sigma[[x]])])}))
            }
            if(length(starting)!=(p*npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
        }

    }
# Optimizer
if(optimization!="subplex"){
estim<-optim(par=starting,fn=function (par) { lik.Mult(par=par, dat=data, C=C2, D=D, index.mat=index.mat,sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint)$loglik },control=control,hessian=TRUE,method=optimization)   #mettre les options maxit et method dans le menu
}else{
estim<-subplex(par=starting,fn=function (par){lik.Mult(par=par,dat=data,C=C2, D=D, index.mat=index.mat, sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint)$loglik},control=control,hessian=TRUE)   #mettre les options maxit et method dans le menu
}

}else if(model=="BM1"){
##---------------------Optimization BM1---------------------------------------##
    if(constraint==FALSE){
        # number of parameters
        npar=(k*(k+1)/2)
        # initial values for the optimizer
        if(is.null(param[["sigma"]])==TRUE){
            sig1<-varBM(tree,data,n,k)
            starting<-sym.unpar(sig1)
        }else{
            if(length(param$sigma)==npar){
                starting<-param$sigma
            }else{
                starting<-sym.unpar(param$sigma)
            }
            if(length(starting)!=(npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
        }
    }else if(constraint=="diagonal"){
##---------------------Optimization BM1 diagonal------------------------------##
npar=k
# initial values for the optimizer
sig1<-varBM(tree,data,n,k)
# starting values following Adams (2012)
sigma.mn<-diag(sig1)
if(is.null(param[["sigma"]])==TRUE){
    starting=c(sigma.mn)
    
}else{
    if(length(param$sigma)==npar){
        starting<-param$sigma
    }else{
        starting<-c(param$sigma[[1]][[1]],param$sigma[[1]][lower.tri(param$sigma[[1]])])
    }
    if(length(starting)!=(npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
    
}

    }else{
##---------------------Optimization BM1 constrained---------------------------##
        #Same as Adams 2012
        npar=(k*(k-1)/2)+1
        # initial values for the optimizer
        sig1<-varBM(tree,data,n,k)
        # starting values following Adams (2012)
        sigma.mn<-mean(diag(sig1))
        R.offd<-rep(0,(k*(k-1)/2))
        if(is.null(param[["sigma"]])==TRUE){
            starting=c(sigma.mn,R.offd)
            
        }else{
            if(length(param$sigma)==npar){
                starting<-param$sigma
            }else{
                starting<-c(param$sigma[[1]][[1]],param$sigma[[1]][lower.tri(param$sigma[[1]])])
            }
            if(length(starting)!=(npar)){stop("The number of starting values for the rates matrix do not match, see ?mvBM for providing user specified starting values")}
            
        }
    }

# Optimizer
if(optimization!="subplex"){
estim<-optim(par=starting,fn=function(par){lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,k=k,constraint=constraint)$loglik},control=control,hessian=TRUE,method=optimization)
}else{
estim<-subplex(par=starting,fn=function(par){lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,k=k,constraint=constraint)$loglik},control=control,hessian=TRUE)
} 
}

##-----------------Summarizing results----------------------------------------##
if(model=="BMM"){
matResults<-buildSigma(estim$par,index.mat,sig,model,constraint)
resultList<-array(dim = c(k, k, p))
states=vector()
 for(i in 1:p){
 resultList[,,i]<-matResults[[i]]
 states[i]<-colnames(tree$mapped.edge)[i]#multi.tre[[i]]$state
}
dimnames(resultList)<-list(colnames(data), colnames(data), states)

#ancestral states estimates
anc<-lik.Mult(par=estim$par,dat=data,C=C2,D=D,index.mat=index.mat,sig=sig,error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method,constraint=constraint)$ancstate
if(param$smean==TRUE){
    anc<-matrix(anc,nrow=1)
    colnames(anc)<-colnames(data)
    rownames(anc)<-"theta"
}else{
    anc<-matrix(anc,nrow=p)
    colnames(anc)<-colnames(data)
    rownames(anc)<-colnames(tree$mapped.edge)
}

}else if(model=="BM1"){
resultList<-buildSigma(estim$par,NULL,NULL,model,constraint)
 colnames(resultList)<-colnames(data)
 rownames(resultList)<-colnames(data)
 #ancestral states estimates
anc<-matrix(lik.BM1(par=estim$par,dat=data,C=C1,D=D,error=error, k=k, n=n, precalcMat=precalcMat, method=method,constraint=constraint)$ancstate,nrow=1)
colnames(anc)<-colnames(data)
rownames(anc)<-"theta"
}

# LogLikelihood
LL<--estim$value
# models parameters
if(model=="BMM"){
    if(param$smean==TRUE){
        nparam=k+length(unique(estim$par))#k+(p*k) = p for each regimes, k for each rates, k for each ancestral states   or(k+length(unique(index.mat))?
    }else{
        nparam=k*p+length(unique(estim$par))
    }
}else if(model=="BM1"){
nparam=k+length(estim$par)        #k+k= k for each rates and k for each ancestral states
}
# AIC
AIC<--2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1989
# Maybe n need to be changed by length(data)? Moreover it can change when there is missing cases
##---------------------Diagnostics--------------------------------------------##

if(estim$convergence==0 & diagnostic==TRUE){  
cat("successful convergence of the optimizer","\n")
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
if(constraint==TRUE | constraint=="diagonal"){
cat("-- Summary results for multiple constrained rates",model,"model --","\n")
}else if(constraint=="correlation"){
cat("-- Summary results for common correlation ",model,"model --","\n")
}else if(constraint=="shared"){
cat("-- Summary results for shared eigenvectors ",model,"model --","\n")
}else if(constraint=="proportional"){
cat("-- Summary results for proportional rates matrices ",model,"model --","\n")
}else{
cat("-- Summary results for multiple rates",model,"model --","\n")
}
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters","\n")
cat("\n")
cat("Estimated rates matrix","\n")
cat("______________________","\n")
print(resultList)
cat("\n")
cat("Estimated root state","\n")
cat("______________________","\n")
print(anc)
cat("\n")
    }
##-------------------Save infos in parameters---------------------------------##
param$model<-model
param$constraint<-constraint
param$nparam<-nparam
param$nbspecies<-n
param$ntraits<-k
param$nregimes<-p
param$method<-method
param$optimization<-optimization
param$traits<-colnames(data)
##-------------------Store results--------------------------------------------##

 results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=anc, sigma=resultList ,convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.value, param=param)


class(results)<-c("mvmorph","bm")
invisible(results)
#End
}
