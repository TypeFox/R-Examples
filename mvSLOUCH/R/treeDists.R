.calcTreeDist<-function(phyltree,termNames){
## the tree is assumed to be in OUCH format
## terminalNames has to be given so that we have correspondece to traits
## function returns a list with three components
## [[1]], $vSpeciesTime : vector of times of current species
## [[2]], $mTreeDist : matrix where M[i,j]= distance from node i to the last common ancestor of i and j
## [[3]], $mAncestorTimes : matrix where M[i,j]= time of most recent common ancestor of nodes i and j
## [[4]], NULL placeholder
## [[5]], solve([[3]])
## Note that this will NOT be a symmetric matric unless the tree is ultrametric (all present nodes on same level)
## this function overwrites all of the previous recursive functions, which are not needed anymore

    lDists<-vector("list",5)
    names(lDists)<-c("vSpeciesTime","mTreeDist","mAncestorTimes","vSpeciesPairs","invmAncestorTimes")


    nterm<-length(termNames)
    ## first entry ------------------------------------------------------------------------------------
    lDists[[1]]<-rep(0,nterm)
    names(lDists[[1]])<-termNames
    for (i in 1:length(termNames)){## calculate each distance
        iNodeNum<- which(phyltree@nodelabels==termNames[i])[1]
        lDists[[1]][termNames[i]]<- phyltree@times[iNodeNum]
    }
    ##-------------------------------------------------------------------------------------------------
    
    ## second entry -----------------------------------------------------------------------------------
    lDists[[2]]<-matrix(0,ncol=nterm,nrow=nterm)
    colnames(lDists[[2]])<-termNames
    rownames(lDists[[2]])<-termNames
    for (i in 1:length(termNames)){
        lDists[[2]][i,i]<-0 ## same will be 0
        iIndex<- which(phyltree@nodelabels==termNames[i])[1]
        if ((i+1)<=length(termNames)){
            for (j in (i+1):length(termNames)){
                jIndex<- which(phyltree@nodelabels==termNames[j])[1]
                vCommonAncestors<-intersect(phyltree@lineages[[iIndex]],phyltree@lineages[[jIndex]])
                LastIndex<-vCommonAncestors[1]
                LastDist<-phyltree@times[LastIndex]
                if (length(vCommonAncestors)>1){
	    	    for (k in 2:length(vCommonAncestors)){## Look for the Last Common Ancestor
			CurrIndex<-vCommonAncestors[k]
			CurrDist<-phyltree@times[CurrIndex]
			if (CurrDist>LastDist){
			    LastIndex<-CurrIndex
			    LastDist<-CurrDist
			}
		    }
		}
                lDists[[2]][i,j]<-phyltree@times[iIndex]-phyltree@times[LastIndex]
                lDists[[2]][j,i]<-phyltree@times[jIndex]-phyltree@times[LastIndex]
    	    }
	}
    }
    ##-------------------------------------------------------------------------------------------------
    lDists[[3]]<-.calcAncestorTimes(lDists[[2]],lDists[[1]])
    lDists[[4]]<-NULL
    lDists[[5]]<-solve(lDists[[3]])
    lDists
}

.calcAncestorTimes<-function(mTreeDist,vSpecDist){
    n<-length(vSpecDist)
    mAncestorTimes<-matrix(1:n^2,nrow=n,ncol=n) ## need to calculate all the times of ancestors from mTreeDist,mSpecDist
    mAncestorTimes<-apply(mAncestorTimes,c(1,2),function(ij,mTreeDist,vSpecDist,n){
                i<-(ij-1)%%n+1
                j<-(ij-1)%/%n+1
                vSpecDist[i]-mTreeDist[i,j]
    },"mTreeDist"=mTreeDist,"vSpecDist"=vSpecDist,"n"=n)
    colnames(mAncestorTimes)<-colnames(mTreeDist)
    rownames(mAncestorTimes)<-rownames(mTreeDist)    
    mAncestorTimes
}

.calculate.Tree.dists<-function(PhylTree,UserTermLabels=NULL){
    TreeTermLabels<-PhylTree@nodelabels[PhylTree@term]
    n<-length(TreeTermLabels)
    lDists<-.calcTreeDist(PhylTree,TreeTermLabels)
    lPrecalc<-vector("list",6)
    names(lPrecalc)<-c("mSpecDist","mTreeDist","mAncestorTimes","vSpeciesPairs","invmAncestorTimes","tree.height")
    ## the presence of mSpecDist is used in .Params.summary to determine if CIs are to be calculated

    lPrecalc$tree.height<-PhylTree@depth
    
    if (is.null(UserTermLabels)){lPrecalc$mSpecDist<-lDists[[1]][TreeTermLabels]}
    else{lPrecalc$mSpecDist<-lDists[[1]][UserTermLabels]}
    if(!(is.matrix(lPrecalc$mSpecDist))){lPrecalc$mSpecDist<-matrix(lPrecalc$mSpecDist,nrow=1)}
    if (is.null(UserTermLabels)){colnames(lPrecalc$mSpecDist)<-TreeTermLabels}
    else{colnames(lPrecalc$mSpecDist)<-UserTermLabels}
    
    if (is.null(UserTermLabels)){lPrecalc$mTreeDist<-lDists[[2]][TreeTermLabels,TreeTermLabels]}
    else{lPrecalc$mTreeDist<-lDists[[2]][UserTermLabels,UserTermLabels]}
    if (is.null(UserTermLabels)){lPrecalc$mAncestorTimes<-lDists[[3]][TreeTermLabels,TreeTermLabels]}
    else{lPrecalc$mAncestorTimes<-lDists[[3]][UserTermLabels,UserTermLabels]}
    lPrecalc$vSpeciesPairs<-matrix(1:n^2,n,n)[upper.tri(matrix(1,n,n),diag=TRUE)]
    if (is.null(UserTermLabels)){lPrecalc$invmAncestorTimes<-lDists[[5]][TreeTermLabels,TreeTermLabels]}
    else{lPrecalc$invmAncestorTimes<-lDists[[5]][UserTermLabels,UserTermLabels]}

    lPrecalc
}
