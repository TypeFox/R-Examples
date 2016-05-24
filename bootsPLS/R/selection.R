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


# random.subsampling(Y)
# prediction.formula(data,H,uloadings,vloadings,CH,means.X,means.Y,sigma.X,sigma.Y)
# pre.screening(X,coeff)

# random subsamplings
component.selection=function(object,alpha)
{
    # distance: "max.dist", "centroids.dist", "mahalanobis.dist"

    #object is a bootsPLS object
    if(missing(object)) stop("`object' is missing")
    if(missing(alpha)) alpha=0.01
    
    if(!any(class(object)=="bootsPLS")) stop("problem class of object")


    #summary of the classification accuracy over all the replications
    MC=apply(object$ClassifResult,3,function(y){apply(y,3,function(x){sum(x)-sum(diag(x))})})
    
    
    max=ncol(object$nbr.var) #number max of components included
    pval=NULL
    opt=1 #initialise the first optimal number of genes
    for(j in 2:max)
    {
        
        pval[j]=t.test(MC[,opt],MC[,j],alternative="greater")$p.value #t.test of "is adding X genes improves the overall results"
        
        if( (pval[j]< (alpha)))
        {
            opt=j #if the p-value is lower than 0.05, the optimal number of genes is updated
            #cat("opt.temp ", opt,"\n") #print the temporary optimal number of genes for MC
        }
    }
    #if(opt.MC==1)
    #cat("opt.temp ", opt,"\n") #print the temporary optimal number of genes for MC

    cat("Number of chosen components:",opt,"\n")

    out=list(pval=pval,opt=opt,object=object,alpha=alpha)
    structure(out,class="component.selection")

}


variable.selection=function(object,ncomp,alpha,limit)
{
    #object:bootsPLS object
    #ncomp: number of optimal component
    #alpha: level of the tests
    #limit: vector of maximal number of genes to include on each component
    # distance: "max.dist", "centroids.dist", "mahalanobis.dist"
    
    if(missing(object)) stop("`object' is missing")
    if(missing(alpha)) alpha=0.01
    if(missing(ncomp)) ncomp=ncol(object$nbr.var)
    if(missing(limit)) limit=rep(ncol(object$frequency),ncomp)-1
    method=object$data$method
    
    
    
    if(!any(class(object)=="bootsPLS")) stop("problem class of object")


    #loop on the number of optimal component
    subsamplings=list()
    pval=keepX.constraint=list()
    opt=NULL
    for(compi in 1:ncomp)
    {
        
        
        variables=sort(object$frequency[compi,],decreasing=TRUE)
        dup=duplicated(variables[limit[compi]:1])[limit[compi]:1]

        first=which(!dup)[1] #give the first subset of genes. e.g 4genes have the same stability 1.00, first=4
        subsamplings[[compi]]=matrix(NA,nrow=100,ncol=limit[compi])
        for(j in (first):limit[compi])
        {
            if(!dup[j])
            {

                out=CI.prediction(X=object$data$X,Y=object$data$Y,keepX.constraint=list(names(variables)[1:j]),ncomp=1)
                
                # match distance and output in out
                ind.method=which(dimnames(out$ClassifResult)[[5]]==method)
                subsamplings[[compi]][,j]=apply(out$ClassifResult[,,,,ind.method],3,function(x){sum(x)-sum(diag(x))})
            }
        }
        colnames(subsamplings[[compi]])=names(variables)[1:limit[compi]]
        
        
        pval[[compi]]=numeric()
        opt.temp=first #initialise the first optimal number of genes
        #cat("comp.",compi," opt.temp ", opt.temp.MC,"\n") #print the temporary optimal number of genes for MC

        for(j in (first+1):limit[compi])
        {
            if(!dup[j])
            {
                pval[[compi]][j]=t.test(subsamplings[[compi]][,opt.temp],subsamplings[[compi]][,j],alternative="greater")$p.value #t.test of "is adding X genes improves the overall results"
                if( pval[[compi]][j]< (alpha) )
                {
                    opt.temp=j #if the p-value is lower than alpha, the optimal number of genes is updated
                    #cat("comp.",compi," opt.temp.MC ", opt.temp.MC,"\n") #print the temporary optimal number of genes for MC
                }
            }
        }
        opt[compi]=opt.temp

    keepX.constraint[[compi]]=names(variables)[1:opt.temp]

    cat("comp.",compi,"- Number of chosen variables:",opt.temp,"\n")

    }
    names(opt)=paste("comp.",1:ncomp,sep="")
    names(pval)=paste("comp.",1:ncomp,sep="")
    names(subsamplings)=paste("comp.",1:ncomp,sep="")
    names(keepX.constraint)=paste("comp.",1:ncomp,sep="")

    out=list(pval=pval,opt=opt,keepX.constraint=keepX.constraint,subsamplings=subsamplings,object=object,alpha=alpha)
    structure(out,class="variable.selection")
    
}
