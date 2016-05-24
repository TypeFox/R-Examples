# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 13-05-2014
# last modification: 15-03-2015
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# perform the replication of the splsda on random subsamplings of the data

bootsPLS=function(X,
        Y,
        near.zero.var,
        many,
        ncomp,
        method = c("max.dist", "centroids.dist", "mahalanobis.dist"),
        save.file,
        ratio,
        kCV,
        grid,
        showProgress=TRUE)
{
    #-------------------------- bootstrapped sPLS-DA and Cross-Validation
    # X input data
    # Y factor input
    # near.zero.var: discard variables with near to zero variance, see function from mixOmics package
    # many= number of resampling. On each subsampling is performed the bootstrap sPLS-DA (with CV to choose the parameters)
    # ncomp= number of components
    # method= which distance should be used to classify the samples?
    # save.file= name of the file with the results to save
    # ratio: proportion of samples left out in the first random subsampling
    # kCV: number of k-fold in the cross validation
    # grid: CV grid, a list of length ncomp. Inside is a list of value of keepX to test in the CV
    # showProgress=TRUE, show the progress of the iteration
    
    if(missing(X)) stop("missing X")
    if(missing(Y)) stop("missing Y")
    if(missing(near.zero.var)) near.zero.var=TRUE
    if(missing(many)) many=50
    if(missing(ncomp)) ncomp=2
    if(missing(kCV)) kCV=10

    check=Check.entry.bootsPLS(X,Y)
    X=check$X
    Y=check$Y

    nlevelY=nlevels(Y)
    
    if(missing(grid))
    {
        grid=list()
        for(i in 1:ncomp)
        grid[[i]]=1:min(40,ncol(X))
    }
    
    if(missing(method)) method="max.dist"
    
    #construct a dummy matrix
    Y.mat=matrix(0,nrow=nrow(X),ncol=nlevelY)
    for(i in 1:nlevelY)
    {
        Y.mat[which(Y==levels(Y)[i]),i]=1
    }
    colnames(Y.mat)=levels(Y)
    
    #remove genes with nearzerovar function
    nzv=list()
    if(near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            names.remove=colnames(X)[nzv$Position]
            warning("Zero- or near-zero variance predictors have been discarded.\n See $nzv for problematic predictors.")
            X = X[, -nzv$Position,drop=FALSE]
            if(ncol(X)==0) {stop("No more predictors")}
        }else{nzv$Position=NULL}
    }
    
    #initialise some parameters that we want to record
    P=ncol(X)
    ClassifResult=array(0,c(nlevelY,nlevelY,ncomp,many))
    rownames(ClassifResult)=levels(Y)
    colnames(ClassifResult)=paste("predicted.as.",levels(Y),sep="")
    dimnames(ClassifResult)[[3]]=paste("comp.",1:ncomp,sep="")
    dimnames(ClassifResult)[[4]]=paste("iteration.",1:many,sep="")
    
    selection.variable=array(0,c(ncomp,P,many))
    dimnames(selection.variable)[[2]]=colnames(X)
    loadings.X=array(0,c(P,ncomp, many))
    dimnames(loadings.X)[[1]]=colnames(X)
    nbr.var=NULL
    
    #--------------------------------------------------------------------------------------------------------
    #---------------------		LOOP   	---------------------
    #--------------------------------------------------------------------------------------------------------
    learning.sample=matrix(0,nrow=nrow(X),ncol=many) #record which sample are in the learning set
    prediction=array(0,c(nrow(X),many,ncomp)) #record the class associated to each sample (either in learning or test set)
    
    for(abc in 1:many)
    {
        if(showProgress)
        cat("iteration ",abc,"\n")
        
        #--------------- 1st step: random subsampling
        
        A=suppressWarnings(random.subsampling(Y,ratio))
        #A contains the sample we want to keep in the learning set, -A in the test set
        learning.sample[A,abc]=1
        
        data.learn.signature=X[A,]
        Y.learn.signature=Y[A]
        Y.mat.learn.signature=Y.mat[A,]
        
        data.test.signature=X[-A,]
        Y.test.signature=Y[-A]
        Y.mat.test.signature=Y.mat[-A,]
        
        #we scale data and Y, the CV is done with the scaled matrices, the prediction will be done on the unscaled matrices using the means and variances
        data.learn.scale=scale(data.learn.signature)
        Y.mat.learn.scale=scale(Y.mat.learn.signature)
        means.X=attr(data.learn.scale,"scaled:center")
        sigma.X=attr(data.learn.scale,"scaled:scale")
        means.Y=attr(Y.mat.learn.scale,"scaled:center")
        sigma.Y=attr(Y.mat.learn.scale,"scaled:scale")
        
        remove=which(sigma.X==0)
        if(length(remove)>0)
        {
            data.learn.scale=data.learn.scale[,-remove]
            data.test.signature=data.test.signature[,-remove]
            data.learn.signature=data.learn.signature[,-remove]
            means.X=means.X[-remove]
            sigma.X=sigma.X[-remove]
        }
        
        nzv.temp=nearZeroVar(data.learn.scale)
        if(length(nzv.temp$Position)>0)
        {
            data.learn.scale=data.learn.scale[,-nzv.temp$Position]
            data.test.signature=data.test.signature[,-nzv.temp$Position]
            data.learn.signature=data.learn.signature[,-nzv.temp$Position]
            means.X=means.X[-nzv.temp$Position]
            sigma.X=sigma.X[-nzv.temp$Position]
        }
        
        #--------------- kCV on each component to find the optimal number of variables per component
        SIGN=matrix(0,nrow=ncomp,ncol=ncol(X))
        CH=NULL
        uloadings.X=NULL # the one of size p = ncol(X)
        ind.var=NULL
        keepX.constraint=list()
        
        stop=0
        for(num.comp in 1:ncomp)
        {
            if(showProgress)
            cat("iteration",abc,"- comp",num.comp,"\n")
            
            text=paste("iteration",abc,"comp",num.comp,"\n")
            if(TRUE) #debug: get rid of the CV part that is time-consuming
            {
                #perform a Cross validation to tune the number of variables to keep on component num.comp
                CV=CV.spls.hybrid(data=data.learn.signature,Y.factor=Y.learn.signature,Y.mat=Y.mat.learn.signature,
                k=kCV,num.comp=num.comp,keepX.constraint=keepX.constraint,grid.comp=grid[[num.comp]],
                method=method,showProgress=showProgress) #gives number of variables on component num.comp and the threshold on the same component
                # data: input data matrix
                # Y.factor: factor Y
                # Y.mat: dummy matrix constructed from Y.factor
                # k: number of folds for the CV
                # num.comp: which component we are tuning
                # keepX.constraint: which variables are kept on the first num.comp-1 components
                # grid.comp: grid of keepX that is to be tested in the CV
                # method= which distance should be used to classify the samples?
                # showProgress=TRUE, show the progress of the iteration
                
                
                nbr.var.opt=CV$nbr.var.opt
                #if(nbr.var.opt==0){save(list=ls(),file=save.file)}
            }else{
                
                ind=sample(1:ncol(data.learn.scale),3)
                nbr.var.opt=length(ind)
            }
            #----------- learning the model with nbr.var.opt
            res=spls.hybrid(data.learn.scale,Y.mat.learn.scale,keepX.constraint=keepX.constraint,ncomp=num.comp,keepX=nbr.var.opt,near.zero.var=FALSE)
            ind=which(res$loadings$X[,num.comp]!=0)
            names(ind)=colnames(data.learn.scale)[ind]
            
            keepX.constraint[[num.comp]]=names(ind)
            
            #record the signature
            signature.value.X=matrix(0,nrow=ncol(X),ncol=1)
            a=match(names(ind),colnames(X))
            signature.value.X[a]=res$loadings$X[ind,num.comp]
            uloadings.X=cbind(uloadings.X,signature.value.X)
            SIGN[num.comp,a]=1
            
            ind.var=c(ind.var,nbr.var.opt)
            
        }# end num.comp
        
        names(keepX.constraint)=paste("comp.",1:ncomp,sep="")
        if(showProgress) {print(keepX.constraint)}


        #learning the model on num.comp component with spls.constraint
        res=spls.hybrid(data.learn.scale,Y.mat.learn.scale,keepX.constraint=keepX.constraint,ncomp=num.comp)
        uloadings=res$loadings$X
        tvariates=res$variates$X
        CH=res$mat.c
        
        
        #----------- test of the signature on the unscaled-learning.set
        out=prediction.formula(X.test=data.learn.signature,ncomp=ncomp,Y.scaled=Y.mat.learn.scale,unmap.Y=Y.mat.learn.signature,
        variates.X=tvariates,uloadings=uloadings,CH=CH,means.X=means.X,means.Y=means.Y,sigma.X=sigma.X,sigma.Y=sigma.Y,method=method)
        
        if(method == "max.dist"){predicted.learn=out$class$max.dist}
        if(method == "centroids.dist"){predicted.learn=out$class$centroids.dist}
        if(method == "mahalanobis.dist"){predicted.learn=out$class$mahalanobis.dist}
        prediction[which(learning.sample[,abc]==1),abc,]=predicted.learn    # record the prediction
        
        #----------- test of the signature on the unscaled-test.set
        out=prediction.formula(X.test=data.test.signature,ncomp=ncomp,Y.scaled=Y.mat.learn.scale,unmap.Y=Y.mat.learn.signature,
        variates.X=tvariates,uloadings=uloadings,CH=CH,means.X=means.X,means.Y=means.Y,sigma.X=sigma.X,sigma.Y=sigma.Y,method=method)
        
        if(method == "max.dist"){predicted=out$class$max.dist}
        if(method == "centroids.dist"){predicted=out$class$centroids.dist}
        if(method == "mahalanobis.dist"){predicted=out$class$mahalanobis.dist}
        prediction[which(learning.sample[,abc]==0),abc,]=predicted    # record the prediction
        
        
        
        #--------record of the classification accuracy for each level of Y
        for(i in 1:nlevelY)
        {
            ind.i=which(Y.test.signature==levels(Y)[i])
            for(j in 1:nlevelY)
            {
                ClassifResult[i,j,,abc]=apply(predicted,2,function(x){sum(x[ind.i]==j)})
                
            }
        }
        
        selection.variable[,,abc]=SIGN
        loadings.X[,,abc]=uloadings.X
        nbr.var=rbind(nbr.var,ind.var)
        
        #calculation of the frequency of selection for each variable, on each components, after the abc replications.
        # this is done to ensured that the file saved can be re-used.
        frequency=matrix(0,nrow=ncomp,ncol=dim(loadings.X)[1])
        for(j in 1:abc)
        {
            for(k in 1:ncomp)
            {a=which(loadings.X[,k,j]!=0)
                frequency[k,a]=frequency[k,a]+1 #add 1 everytime the gene is selected
            }
        }
        frequency=frequency/abc #get the probability of selection (percentage of times each gene is selected, per component
        colnames(frequency)=colnames(X)
        
        out=list(ClassifResult=ClassifResult,loadings.X=loadings.X,selection.variable=selection.variable,frequency=frequency,
        nbr.var=nbr.var,learning.sample=learning.sample,prediction=prediction,data=list(X=X,Y=Y,method=method),nzv=nzv)
        structure(out,class="bootsPLS")

        data=list(X=X,Y=Y,method=method)
        if(!missing(save.file))
        save(ClassifResult,loadings.X,selection.variable,frequency,nbr.var,learning.sample,prediction,data,nzv,file=save.file)
        
    }
    out
    structure(out,class="bootsPLS")

}#end function

