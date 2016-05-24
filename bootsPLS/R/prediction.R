# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 26-05-2014
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


# perform the prediction on a learning set and a testing dataset


prediction=function(object,X,Y,keepX.constraint,ncomp,X.test,
            CI,many,subsampling.matrix,ratio,level.CI,
            save.file)
{
    #object: bootsPLS.model object, from fit.model. If `object' is given, only X.test is used
    
    # if object is missing, X, Y, keepX.constraint and X.test are used
    
    # X=learning data, keepX.constraint must be included in colnames(X)
    # Y=learning Y
    # keepX.constraint= the signature, a list of genes of length ncomp or their position in the matrix X
    # ncomp= number of components
    # X.test
    
    # CI: logical, should Confidence Interval be calculated?
    # many: how many subsampling for the CI
    # subsampling.matrix: Gives the samples to subsample as an internal learning set
    # ratio: Number between 0 and 1. It is the proportion of the n samples that are put aside and considered as an internal testing set. The (1-ratio)*n samples are used as a training set
    # level.CI: A 1- level.CI% confidence interval is calculated
    
    
    if(missing(object))
    {
        # check and fit the model with X, Y, keepX.constraint, ncomp
        if(missing(X)) stop("missing X")
        if(missing(Y)) stop("missing Y")
        if(missing(keepX.constraint)) stop("missing keepX.constraint")
        if(missing(ncomp)) ncomp=length(keepX.constraint)
        if(ncomp>length(keepX.constraint)) stop("The number of components has to be lower than the length of keepX.constraint")
        if(ncomp!=length(keepX.constraint)) {keepX.constraint.temp=keepX.constraint; for(i in (ncomp+1):length(keepX.constraint)) keepX.constraint.temp[[i]]=NULL ; keepX.constraint=keepX.constraint.temp}
   
       
        #fit the model
        object=fit.model(X=X,Y=Y,ncomp=ncomp,keepX.constraint=keepX.constraint)# spls.constraint object
    }else{
        if(!any(class(object)=="spls.constraint")) stop("problem class of object")
    }
    
    if(missing(CI)) CI=FALSE
    
    #use the parameters from the model
    
    data=object$X
    Y=object$Y #factor
    keepX.constraint=object$data$keepX.constraint
    ncomp=length(keepX.constraint)
    
    #data=X[,unique(unlist(keepX.constraint)),drop=FALSE] #pick the genes
    
    Y.mat=object$data$Y.mat
    Y.mat.scale=object$Y

    if(length(object$nzv$Position)>0)
    {
        means.X=object$coeff$means.X[-object$nzv$Position]
        sigma.X=object$coeff$sigma.X[-object$nzv$Position]
        
    }else{
        means.X=object$coeff$means.X
        sigma.X=object$coeff$sigma.X
    }
    
    means.Y=object$coeff$means.Y
    sigma.Y=object$coeff$sigma.Y

        
    
    uloadings=object$loadings$X
    tvariates=object$variates$X
    CH=object$mat.c
    
    #----------- test of the signature on the unscaled-training.set
    out=prediction.formula(X.test=data,ncomp=ncomp,Y.scaled=Y.mat.scale,unmap.Y=Y.mat,
    variates.X=tvariates,uloadings=uloadings,CH=CH,means.X=means.X,means.Y=means.Y,sigma.X=sigma.X,sigma.Y=sigma.Y)
    predicted.learn=out$class
    Y.hat.learn=out$Y.hat
    
    out=list(object=object,predicted.learn=predicted.learn,Y.hat.learn=Y.hat.learn)#,t.learn=t.learn)


    #----------- test of the signature on the unscaled-test.set
    if(!missing(X.test))
    {
        #check entry
        check=Check.entry.X(X.test)
        X.test.temp=check$X
        
        # prep X.test
        match.indice=match(colnames(data),colnames(X.test.temp))
        missing.genes=sum(is.na(match.indice))
        if(sum(is.na(match.indice))>0)
        {
            warning("Some genes are missing from the signature")
        }
        if(sum(is.na(match.indice))==ncol(data))
        {
            warning("All genes are missing from the signature, prediction shouldn't be trusted")
        }
        
        X.test.temp=X.test.temp[,match.indice,drop=FALSE]
        X.test.temp[is.na(X.test.temp)]=0


        # make the prediction
        out.temp=prediction.formula(X.test=X.test.temp,ncomp=ncomp,Y.scaled=Y.mat.scale,unmap.Y=Y.mat,
        variates.X=tvariates,uloadings=uloadings,CH=CH,means.X=means.X,means.Y=means.Y,sigma.X=sigma.X,sigma.Y=sigma.Y)
        
        #record the result
        predicted.test=out.temp$class
        Y.hat.test=out.temp$Y.hat
        #t.pred=out$t.hat
        

         
         out$X.test=X.test.temp
         out$predicted.test=predicted.test
         out$Y.hat.test=Y.hat.test
         out$missing.genes=missing.genes
         
         if(CI==TRUE)
         {
             out.CI=CI.prediction(object=object,keepX.constraint=keepX.constraint,ncomp=ncomp,many=many,subsampling.matrix=subsampling.matrix,
             ratio=ratio,X.test=X.test,level.CI=level.CI,save.file=save.file)
             out$out.CI=out.CI
  
             if(!missing(save.file)) save(object,X.test,predicted.learn,Y.hat.learn,#t.learn,
             Y.hat.test,predicted.test,out.CI,#t.pred,
             file=save.file)
         }else{
             
             if(!missing(save.file)) save(object,X.test,predicted.learn,Y.hat.learn,#t.learn,
             Y.hat.test,predicted.test,#t.pred,
             file=save.file)
         }



    }else{
        if(!missing(save.file)) save(object,predicted.learn,#t.learn,
            Y.hat.learn,file=save.file)
    }
    
    out
}#end function









