################################################################################
##                                                                            ##
##                               mvMORPH: mv.precalc                          ##
##                                                                            ##
##  Created by Julien Clavel - 01-11-2014                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex                                 ##
##                                                                            ##
################################################################################

mv.Precalc<-function(tree, nb.traits=1, scale.height=FALSE, param=list(pivot="MMD", method=c("sparse"), smean=TRUE, model="OUM")){
    
    if(inherits(tree,"phylo") & inherits(tree,"multiPhylo")){ ## A MODIFIER???
        stop("The tree is not a \"phylo\" object")
    }

    n<-length(tree$tip.label)
    #scale the tree
    if(scale.height==TRUE){
        maxHeight<-max(nodeHeights(tree))
        tree$edge.length<-tree$edge.length/maxHeight
        tree$mapped.edge<-tree$mapped.edge/maxHeight
    }
    #set data as a matrix if a vector is provided instead
    #if(!is.matrix(data)){data<-as.matrix(data)}
    #if(is.vector(data)){k<-1 }else{ k<-ncol(data)}
    p<-nb.traits
    #set tree order
    ind<-reorder.phylo(tree,order="postorder", index.only=TRUE)
    tree$edge<-tree$edge[ind,]
    tree$edge.length<-tree$edge.length[ind]


    #Is SIMMAP?
    if(!is.null(tree[["mapped.edge"]])){
    tree$mapped.edge<-tree$mapped.edge[ind,]
    tree$maps<-tree$maps[ind]
    }
    
    #Method?
    if(is.null(param[["method"]])){
        param$method="none"
    }
    
    if(param$method!="pic"){
    #compute the variance-covariance matrix
    C1<-vcv.phylo(tree)
    


if(!is.null(tree[["mapped.edge"]])){# & param[["model"]]!="OUM" & param[["model"]]!="OU1"){
    #Build a list of VCV from SIMMAP trees
    multi.tre<-list()
    class(multi.tre)<-"multiPhylo"
    C2<-list()
    for(i in 1:ncol(tree$mapped.edge)){
        multi.tre[[i]]<-tree
        multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
        multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
        temp<-vcv.phylo(multi.tre[[i]])
        # Should we provide the data object?
        # if(any(tree$tip.label==rownames(data))) {
        #    C2[[i]]<-temp[rownames(data),rownames(data)]
        # }else{
            C2[[i]]<-temp
        #}
    }
    }else{
     C2<-NULL
    }
 
    #compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    D<-multD(tree,p,n,smean=param$smean)

    
    ##-----------------------Precalculate regime indexation-----------------------##
    if(!is.null(param[["model"]])){
        
        # Pull (alpha) matrix decomposition
        if(is.null(param[["decomp"]])){
            decomp<-param$decomp<-"symmetricPositive"
        }else{
            decomp<-param$decomp[1]
        }
        # option for computing the variance covariance matrix
        if(is.null(param[["vcv"]])==TRUE){
            if(param$method=="sparse"){
                vcvtype<-"sparse"
            }else{
                if(is.ultrametric(tree)==TRUE){
                    vcvtype<-"ouch"
                        }else{
                    vcvtype<-"mvmorph"
             
                }
            }
            
        }else{
            vcvtype<-param$vcv
        }
        
        if(vcvtype!="sparse" & param$method=="sparse"){
            vcvtype<-"sparse"
            cat("Only \"sparse\" VCV could be used with the \"sparse\" method. See ?mvOU","\n")
            #method<-"sparse"
        }
        
        # root estimation
        if(is.null(param[["root"]])!=TRUE){
            if(param[["root"]]==TRUE || param[["root"]]==FALSE || param[["root"]]=="stationary"){
                root<-param$root
            }else{
                stop("Only TRUE,FALSE or \"stationary\" are accepted for the \"root\" argument in \"param\"")
            }
            if(param[["root"]]=="stationary" & param$model=="OU1"){
                root<-param$root<-FALSE
            }
        }else{
            root<-TRUE
        }
        
        
        if(param$model=="OU1"){
            k<-1
        }else{
            k<-length(colnames(tree$mapped.edge))
        }
    # root to tip lineage indexation
    root2tip <- .Call("seq_root2tipM", tree$edge, n, tree$Nnode)
    # Si OU1 sur un objet 'phylo'
    if(param$model=="OU1"){

                valLineage<-sapply(1:n,function(z){rev(unlist(
                    sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$edge.length[val]<-tree$edge.length[val]},simplify=FALSE)))
                } ,simplify=FALSE)
          
    }else{
        # Donnees de temps par regimes et par branches
    
            valLineage<-sapply(1:n,function(z){rev(unlist(
                sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$maps[[val]]<-tree$maps[[val]]},simplify=FALSE)))
            } ,simplify=FALSE)
    
    }
    
    # Indexer les regimes
    if(param$model=="OUM"){
        if(root==FALSE | root=="stationary"){
            facInd<-factor(colnames(tree$mapped.edge))
        }else if(root==TRUE){
            facInd<-factor(c("_root_state",colnames(tree$mapped.edge)))
        }
        
        indice<-lapply(1:n,function(z){rev(unlist(lapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); factor(names(tree$maps[[val]]),levels=facInd)})))})
        
    }else if(param$model=="OU1"){
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
    
    }else{
    listReg<-NULL
    epochs<-NULL
    }
    }else{
        listReg<-NULL
        epochs<-NULL
        C1<-NULL
        C2<-NULL
        model<-"none"
    }# end if pic
    ##-----------------------Precalculate sparse method-----------------------##
    if(param$method=="sparse"){
        #require(spam)
        # Temporary multivariate VCV
        V<-kronecker((matrix(1,p,p)+diag(p)), C1)
        
        # Check for missing cases? need to provides the data
        #if(any(is.na(data))){
        #        Indice_NA<-which(is.na(as.vector(data)))
        #        V<-V[-Indice_NA,-Indice_NA]
        #        }

        # spam object
        V1<-as.spam(V);
        # precal the cholesky
        if(is.null(param[["pivot"]])){pivot<-"MMD"}else{pivot<-param$pivot}
        ch<-chol(V1,pivot=pivot)
        # Yale Sparse Format indices
        JAr<-V1@colindices-1
        IAr<-V1@rowpointers-1
        ldv<-dim(V)[1]
        resume<-(length(JAr)*100)/(ldv*ldv)
       cat("The density of the matrix is:",resume,"%","\n")
    }else{
        JAr<-NULL
        IAr<-NULL
        ch<-NULL
        V1<-NULL
    }
    ##------------------param results---------------------------------------------##
    param$root<-root
    param$nbtraits<-nb.traits
    ##------------------List results----------------------------------------------##
    
    results<-list(tree=tree,ch=ch, V=V1, JAr=JAr, IAr=IAr, C1=C1, C2=C2, D=D, listReg=listReg, epochs=epochs, model=param$model, param=param)
    class(results)<-"mvmorph.precalc"
    invisible(results)
}
