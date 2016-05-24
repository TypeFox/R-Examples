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


# fit the constraint spls model and record the prediction on many random subsampling with a set of variables (keepX.constraint)


CI.prediction=function(object,X,Y,keepX.constraint,ncomp,many,subsampling.matrix,ratio,X.test,level.CI,save.file)
{
    #object: spls.constraint object, from fit.model. If `object' is given, X,Y,keepX.constraint,ncomp are not used

    # X input
    # Y factor input
    # keepX.constraint= the signature, a list of genes of length ncomp or their position in the matrix X
    # many= #number of resampling. On each subsmapling is performed the bootstrap sPLS-DA (with CV to choose the parameters)
    # save.file= name of the file to save

    if(missing(many)) many=100
    if(missing(level.CI)) level.CI=0.05

    if(missing(object))
    {
        if(missing(X)) stop("missing X")
        if(missing(Y)) stop("missing Y")
        #check the entries
        check=Check.entry.bootsPLS(X,Y)
        X=check$X
        Y=check$Y
        
        if(missing(keepX.constraint)) stop("missing keepX.constraint")
        if(missing(ncomp)) ncomp=length(keepX.constraint)
        if(ncomp>length(keepX.constraint)) stop("The number of components has to be lower than the length of keepX.constraint")
        if(sum(is.na(match(unlist(keepX.constraint),colnames(X))))>0)
        {
            a=sum(is.na(match(unlist(keepX.constraint),colnames(X))))
            stop(paste("keepX.constraint contains",a,"variable(s) that are missing from X"))
        }
        
        if(ncomp!=length(keepX.constraint)) {keepX.constraint.temp=keepX.constraint; for(i in (ncomp+1):length(keepX.constraint)) keepX.constraint.temp[[i]]=NULL ; keepX.constraint=keepX.constraint.temp}
        

        
    }else{
        if(!any(class(object)=="spls.constraint")) stop("problem class of object")
        X=object$data$X
        Y=object$data$Y #factor
        keepX.constraint=object$data$keepX.constraint
        ncomp=length(keepX.constraint)
    }


    X.indice=X[,unique(unlist(keepX.constraint)),drop=FALSE] #pick the genes
    nlevelY=nlevels(Y)

    #construct a dummy matrix
    Y.mat=matrix(0,nrow=nrow(X),ncol=nlevelY)
    for(i in 1:nlevelY)
    {
        Y.mat[which(Y==levels(Y)[i]),i]=1
    }

    colnames(Y.mat)=levels(Y)



    if(missing(X.test))
    {make.prediction=FALSE}else{make.prediction=TRUE}

    #prediction of all the samples
    out=prediction(X=X,Y=Y,keepX.constraint=keepX.constraint,ncomp=ncomp) #this removes the column without variance
    num.distance=length(out$predicted.learn) # gives the number of different distances


    #initialise some parameters that we want to record

    P=ncol(X.indice)
    ClassifResult=array(0,c(nlevelY,nlevelY,ncomp,many,num.distance))
    rownames(ClassifResult)=levels(Y)
    colnames(ClassifResult)=paste("predicted.as.",levels(Y),sep="")
    dimnames(ClassifResult)[[3]]=paste("comp.",1:ncomp,sep="")
    dimnames(ClassifResult)[[4]]=paste("iteration.",1:many,sep="")
    dimnames(ClassifResult)[[5]]=names(out$predicted.learn)
    
    loadings.X=array(0,c(P,ncomp, many))
    rownames(loadings.X)=colnames(X.indice)


    #CH.tot=array(0,c(P,ncomp, many))

    means.X=matrix(0,nrow=many,ncol=P)
    sigma.X=matrix(0,nrow=many,ncol=P)
    means.Y=matrix(0,nrow=many,ncol=nlevelY)
    sigma.Y=matrix(0,nrow=many,ncol=nlevelY)


    learning.sample=matrix(0,nrow=nrow(X.indice),ncol=many) #record which sample are in the learning set
    prediction.X=array(0,c(nrow(X.indice),many,ncomp,num.distance)) #record the class associated to each sample (either in learning or test set)
    dimnames(prediction.X)[[4]]=names(out$predicted.learn)
    if(make.prediction==TRUE)
    {
        prediction.X.test=array(0,c(nrow(X.test),many,ncomp,num.distance)) #record the class associated to each sample (either in learning or test set)
        dimnames(prediction.X.test)[[4]]=names(out$predicted.learn)
        Y.hat.test=array(0,c(nrow(X.test),nlevelY,ncomp,many)) #record the prediction values associated to each sample in test set
    }
    
    for(abc in 1:many)
    {
        #print(paste("iteration",abc))
        
        
        #--------------- 1st step: random subsampling
        if(missing(subsampling.matrix))
        {
            subsampling=random.subsampling(Y,ratio)
        }else{
            subsampling=subsampling.matrix[,abc]
            }
        #subsampling contains the sample we want to keep in the learning set, -subsampling in the test set
        
        learning.sample[subsampling,abc]=1
        
        data.learn.signature=X[subsampling,,drop=FALSE]
        Y.learn.signature=Y[subsampling]
        Y.mat.learn.signature=Y.mat[subsampling,]
        
        data.test.signature=X[-subsampling,,drop=FALSE]
        Y.test.signature=Y[-subsampling]
        
        out=prediction(X=data.learn.signature,Y=Y.learn.signature,keepX.constraint=keepX.constraint,ncomp=ncomp,X.test=data.test.signature) #this removes the column without variance
        
        for(distance in 1:num.distance)
        {
            prediction.X[subsampling,abc,,distance]=out$predicted.learn[[distance]]
            prediction.X[-subsampling,abc,,distance]=out$predicted.test[[distance]]
            
        }
        

        for(distance in 1:num.distance)
        {
            predicted=out$predicted.test[[distance]]
            #--------record of the classification accuracy for each level of Y
            for(i in 1:nlevelY)
            {
                ind.i=which(Y.test.signature==levels(Y)[i])
                for(j in 1:nlevelY)
                {
                    ClassifResult[i,j,,abc,distance]=apply(predicted,2,function(x){sum(x[ind.i]==j)})
                    
                }
            }
            
        }

        #----------- record of the coefficients used in the prediction function
        ind.commun=match(rownames(out$object$loadings$X),rownames(loadings.X))
        loadings.X[ind.commun,,abc]=out$object$loadings$X
        
        if(length(out$object$nzv$Position)>0)
        {
            means.X[abc,ind.commun]=out$object$coeff$means.X[-out$object$nzv$Position]
            sigma.X[abc,ind.commun]=out$object$coeff$sigma.X[-out$object$nzv$Position]
            
        }else{
            means.X[abc,ind.commun]=out$object$coeff$means.X
            sigma.X[abc,ind.commun]=out$object$coeff$sigma.X
        }
        
        
        means.Y[abc,]=out$object$coeff$means.Y
        sigma.Y[abc,]=out$object$coeff$sigma.Y
        
        #----------- record of the prediction value for the test.set
        if(make.prediction==TRUE)
        {
            out=prediction(X=data.learn.signature,Y=Y.learn.signature,keepX.constraint=keepX.constraint,ncomp=ncomp,X.test=X.test) #need to remove somewhere the column without variance
            predicted=out$predicted.test#$max.dist

            for(distance in 1:num.distance)
            {
                prediction.X.test[,abc,,distance]=predicted[[i]]
            }

            
            Y.hat.test[,,,abc]=out$Y.hat.test
            
            if(sum(is.na(out$Y.hat.test))>0) stop("bla")
            
        }



    } #end abc




    if(make.prediction==TRUE)
    {
        dimnames(Y.hat.test)[[1]]=rownames(X.test)
        dimnames(Y.hat.test)[[2]]=levels(Y)
        dimnames(Y.hat.test)[[3]]=paste("comp.",1:ncomp,sep="")
        dimnames(Y.hat.test)[[4]]=paste("iteration.",1:many,sep="")
        
        min.CI.pred=list()
        for(i in 1:nlevelY)
        {
            min.CI.pred[[i]]=apply(Y.hat.test[,i,,,drop=FALSE],c(1,3),function(x){quantile(x,level.CI/2)})
        }
        names(min.CI.pred)=levels(Y)
        
        
        
        max.CI.pred=list()
        for(i in 1:nlevelY)
        {
            max.CI.pred[[i]]=apply(Y.hat.test[,i,,,drop=FALSE],c(1,3),function(x){quantile(x,1-level.CI/2)})
            
        }
        names(max.CI.pred)=levels(Y)

        #CI contains the confidence interval for each sample, per component, per level.
        CI=list()
        for(i in 1:ncomp)
        {
            CI[[i]]=list()
            for(j in 1:nlevelY)
            {
                temp=cbind(min.CI.pred[[j]][,i],max.CI.pred[[j]][,i])
                colnames(temp)=c("lwr","upr")
                CI[[i]][[j]]=temp
                
            }
            names(CI[[i]])=levels(Y)

        }
        names(CI)=paste("comp.",1:ncomp,sep="")

        out=list(CI=CI,Y.hat.test=Y.hat.test,ClassifResult=ClassifResult,
        loadings.X=loadings.X,prediction.X=prediction.X,prediction.X.test=prediction.X.test,learning.sample=learning.sample,
        data=list(X=X,Y=Y,keepX.constraint=keepX.constraint),coeff=list(means.X=means.X,sigma.X=sigma.X,means.Y=means.Y,sigma.Y=sigma.Y))

      }else{
          CI=NULL
          Y.hat.test=NULL
          
          out=list(CI=CI,Y.hat.test=Y.hat.test,ClassifResult=ClassifResult,
          loadings.X=loadings.X,prediction.X=prediction.X,learning.sample=learning.sample,
          data=list(X=X,Y=Y,keepX.constraint=keepX.constraint),coeff=list(means.X=means.X,sigma.X=sigma.X,means.Y=means.Y,sigma.Y=sigma.Y))

    }
      
      if(!missing(save.file)){
          save(prediction.X,learning.sample,loadings.X,
          ClassifResult,X,Y,keepX.constraint,means.X,sigma.X,means.Y,sigma.Y,
          CI,Y.hat.test,
          file=save.file)
      }
      

    out

}
