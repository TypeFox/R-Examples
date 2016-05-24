# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 28-05-2014
# last modification: 10-10-2014
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


# perform the cross validation to find the optimal number of variables to keep on each components


CV.spls.hybrid=function(data,Y.factor,Y.mat,k,num.comp,keepX.constraint,grid.comp,method,showProgress) #gives number of variable on component num.comp
{
    
    # data: input data matrix
    # Y.factor: factor Y
    # Y.mat: dummy matrix constructed from Y.factor
    # k: number of folds for the CV
    # num.comp: which component we are tuning
    # keepX.constraint: which variables are kept on the first num.comp-1 components
    # grid.comp: grid of keepX that is to be tested in the CV
    # method= which distance should be used to classify the samples?
    # showProgress=TRUE, show the progress of the iteration
    
    
    if(missing(grid.comp))
    {
        grid.comp=seq(1,min(40,ncol(data)),1)
    }
    
    
    #-----------k-fold
    
    n=nrow(data)
    nk=n/k
    #---- construction of the k-CV splits while keeping enough samples of every levels in each split / stratified bootstrap
    
    for(i in 1:nlevels(Y.factor))
    {
        ai=sample(which(Y.factor==levels(Y.factor)[i]),replace=FALSE) # random sampling of the samples from level i
        aai=suppressWarnings(split(ai,factor(1:min(k,length(ai)))))                       # split of the samples in k-folds
        if(length(ai)<k)                                                # if one level doesn't have at least k samples, the list is completed with "integer(0)"
        {
            for(j in (length(ai)+1):k)
            aai[[j]]=integer(0)
        }
        assign(paste("aa",i,sep="_"),sample(aai,replace=FALSE))         # the `sample(aai)' is to avoid the first group to have a lot more data than the rest
    }
    
    # combination of the different split aa_i into SAMPLE
    SAMPLE=list()
    for(j in 1:k)
    {
        SAMPLE[[j]]=integer(0)
        for(i in 1:nlevels(Y.factor))
        {
            SAMPLE[[j]]=c(SAMPLE[[j]],get(paste("aa",i,sep="_"))[[j]])
        }
    }# SAMPLE is a list of k splits
    
    
    
    PRED=INDICE=matrix(0,nrow=length(grid.comp),ncol=k)
    
    for (h in 1:k) #loop on the number of part
    {
        if(showProgress==TRUE) {print(h)}
        
        # ------- construction of the learning set
        X.learn.CV=data[-SAMPLE[[h]],]
        Y.learn.CV.mat=Y.mat[-SAMPLE[[h]],]
        
        X.learn.CV.scale=scale(X.learn.CV)
        Y.learn.CV.mat.scale=scale(Y.learn.CV.mat)
        means.X = attr(X.learn.CV.scale, "scaled:center")
        sigma.X = attr(X.learn.CV.scale, "scaled:scale")
        means.Y = attr(Y.learn.CV.mat.scale, "scaled:center")
        sigma.Y = attr(Y.learn.CV.mat.scale, "scaled:scale")
        
        
        TOT=matrix(0,nrow=ncol(X.learn.CV),ncol=length(grid.comp))#record the selected variables
        rownames(TOT)=colnames(X.learn.CV)
        for(ij in 1:length(grid.comp))
        {
            var.test=grid.comp[ij]
            sp.i=try(spls.hybrid(X.learn.CV.scale,Y.learn.CV.mat.scale,keepX.constraint=keepX.constraint,ncomp=num.comp,keepX=var.test,near.zero.var=FALSE),TRUE)
            if(any(class(sp.i)!="try-error"))
            {
                indice=which(sp.i$loadings$X[,num.comp]!=0)
                indice=names(indice)
                a=match(indice,colnames(X.learn.CV))
                TOT[a,ij]=1 #record the selected variables
            }
        }
        
        #------ test on test.set
        X.test.CV=data[SAMPLE[[h]],]
        Y.test.CV=Y.factor[SAMPLE[[h]]]
        
        ClassifResult=array(0,c(nlevels(Y.factor),nlevels(Y.factor)))
        rownames(ClassifResult)=levels(Y.factor)
        colnames(ClassifResult)=paste("predicted.as.",levels(Y.factor),sep="")
        

        pred=indice=matrix(0,length(grid.comp))
        for(j in 1:length(grid.comp))
        {
            #learn the model with a pls on the selected variables, and test it on the test set
            ind=rownames(TOT)[which(TOT[,j]==1)]
            indice[j]=length(ind)
            
            keepX.constraint.CV=keepX.constraint
            keepX.constraint.CV[[num.comp]]=ind
            
            
            #perform a 1dim-spls on the selected variables and deflate the matrices X and Y
            res.CV=try(spls.hybrid(X=X.learn.CV.scale,Y=Y.learn.CV.mat.scale,keepX.constraint=keepX.constraint.CV,ncomp=num.comp),FALSE)
            
            if(FALSE) #debug: save everything to have a look at what's wrong
            {
                if(any(class(res.CV)=="try-error"))
                {
                    save(list=ls(),file="temp.Rdata")
                    stop("error res.CV")
                }
            }
            
            tvariates.temp=res.CV$variates$X
            uloadings.temp=res.CV$loadings$X
            CH.temp=res.CV$mat.c
            
            
            # prediction of X.test.CV (undeflated data matrix) with the CV.parameters, to be compared to Y.test.CV (undeflated data matrix)
            out=prediction.formula(X.test=X.test.CV,ncomp=num.comp,Y.scaled=Y.learn.CV.mat.scale,unmap.Y=Y.learn.CV.mat,
            variates.X=tvariates.temp,uloadings=uloadings.temp,CH=CH.temp,means.X=means.X,means.Y=means.Y,sigma.X=sigma.X,sigma.Y=sigma.Y) #useless calculation of the prediction on the num.comp-1 first components
            
            if(method == "max.dist"){predicted=out$class$max.dist[,num.comp]}
            if(method == "centroids.dist"){predicted=out$class$centroids.dist[,num.comp]}
            if(method == "mahalanobis.dist"){predicted=out$class$mahalanobis.dist[,num.comp]}
            
            
            #--------record of the classification accuracy for each level of Y
            for(i in 1:nlevels(Y.factor))
            {
                ind.i=which(Y.test.CV==levels(Y.test.CV)[i])
                for(ij in 1:nlevels(Y.factor))
                {
                    ClassifResult[i,ij]=sum(predicted[ind.i]==ij)
                    
                }
            }
            
            #calculation of the BER
            ClassifResult.temp=ClassifResult
            diag(ClassifResult.temp)=0
            BER=sum(apply(ClassifResult.temp,1,sum,na.rm=TRUE)/apply(ClassifResult,1,sum,na.rm=TRUE),na.rm=TRUE)/nlevels(Y.factor)
            if(is.na(BER)) BER=0
            
            pred[j]=BER
            
        }#end grid.comp
        PRED[,h]=pred
        INDICE[,h]=indice
    }#fin k
    
    Err_pred=apply(PRED,1,mean,na.rm = TRUE)
    num.var=apply(INDICE,1,mean,na.rm = TRUE)
    ind.mini=which(Err_pred==min(Err_pred),arr.ind=TRUE)

    A=num.var[ind.mini]#gives us the average number of variables used to get the prediction
    # I say: let's look for the minimum obtained with the lesser number of variables
    nbr.genes=which.min(A)
    nbr.var.opt=grid.comp[ind.mini[nbr.genes]]
    return(list(nbr.var.opt=nbr.var.opt,PRED=PRED,INDICE=INDICE,Err_pred=Err_pred,num.var=num.var))
}#end function
