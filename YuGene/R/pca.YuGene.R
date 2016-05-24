# Author : F.Rohart
# created 07-07-2014
# last modified 07-07-2014
#
# pca.YuGene aims at performing a PCA on the data centered by study, using the PCA function from mixOmics

#saving the pca function from mixomics. It's gonna be replaced by a S3 method
pca.mixOmics=pca

#pca is now a S3 method
pca=function(X,...){UseMethod("pca")}

#the default pca is the function from mixOmics
pca.default=function(X,ncomp = 2,
        center = TRUE,
        scale = FALSE,
        max.iter = 500,
        tol = 1e-09,...)
{
    out=pca.mixOmics(X,ncomp,center,scale,max.iter,tol)
}

#PCA for the class YuGene, the difference is the addition of the `study' parameter that allows to center the data per study
pca.YuGene=function(X,study,ncomp = 2,center = TRUE,scale = FALSE,max.iter = 500,tol = 1e-09,...)
{
    if(!missing(study))
    {
        # print("pca.YuGene") # debug:check that the right pca is used
        
        M=length(levels(study))   # number of groups
        
        if(min(table(study))<4) # have to be careful if there's only three samples: centering is not a good idea in that case
        {
            warnings("At least one study has less than 3 samples")
        }
        
        # check that rownames(x) is not NULL
        if(length(rownames(X))==0) rownames(X)=1:nrow(X)
        if(length(unique(rownames(X)))!=nrow(X)) stop('samples should have a unique identifier/rowname')
        
        # split the data
        study.data = study_split(X,study)
        X.list.study = study.data$X.list.study
        
        # center data per study, and concatene the data
        concat.X = NULL
        X.list.study.center = list()
        for (m in 1:M)
        {
            X.list.study.center[[m]] = scale(X.list.study[[m]], center = TRUE,scale=FALSE)
            concat.X = rbind(concat.X, unlist(X.list.study.center[[m]]))
        }
        
        # rename cols of concatenated centered data
        colnames(concat.X) = colnames(X)
        
        # sort the samples as in the original X
        indice.match=match(rownames(X),rownames(concat.X))
        concat.X=concat.X[indice.match,]
        
        # perform the pca from mixOmics on the centered data
        X.center.study=concat.X
        out=pca.default(X=X.center.study,ncomp,center,scale,max.iter,tol)
        
    }else{
        # print("pca.default") # debug:check that the right pca is used

        # perform the pca on the original (non-centered) data
        out=pca.default(X=X,ncomp,center,scale,max.iter,tol)
    }
    
    
    invisible(out)
} # end function
