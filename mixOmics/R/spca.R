# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Fangzhou Yao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia

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



##==========================SPARSE PCA =========================##

spca <- 
function(X, 
         ncomp = 3, 
         center = TRUE, 
         scale = TRUE,
         keepX = rep(ncol(X), ncomp),
         max.iter = 500, 
         tol = 1e-06)
{

    #--scaling the data--#
    X=scale(X,center=center,scale=scale)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance.")

    # check that the user did not enter extra arguments #
    # --------------------------------------------------#
    # what the user has entered
    match.user =names(match.call())
    # what the function is expecting
    match.function = c('X', 'ncomp', 'center', 'scale', 'keepX', 'max.iter', 'tol')
    
    #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
    if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotVar')
    
    
    #--initialization--#
    X=as.matrix(X)
    X.temp=as.matrix(X)
    n=nrow(X)
    p=ncol(X)
    
    # put a names on the rows and columns
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")

    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) X.names = 1:nrow(X)

    if (length(keepX) != ncomp) 
            stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepX > p)) 
            stop("each component of 'keepX' must be lower or equal than ", p, ".")

    if (ncomp > min(ncol(X), nrow(X)))
    stop("use smaller 'ncomp'", call. = FALSE)

    vect.varX=vector(length=ncomp)
    names(vect.varX) = c(1:ncomp)

    vect.iter=vector(length=ncomp)
    names(vect.iter) = c(1:ncomp)

    vect.keepX=vector(length=ncomp)
    names(vect.keepX) = c(1:ncomp)

# KA: to add if biplot function (but to be fixed!)
    #sdev = vector(length = ncomp)

    mat.u=matrix(nrow=n, ncol=ncomp)
    mat.v=matrix(nrow=p, ncol=ncomp)
    colnames(mat.u)=c(1:ncomp)
    colnames(mat.v)=c(1:ncomp)
    rownames(mat.v)=colnames(X)

    #--loop on h--#
    for(h in 1:ncomp){
       
       #--computing the SVD--#
       svd.X=svd(X.temp)
       u.new = svd.X$u[,1] 
       v.new = svd.X$d[1]*svd.X$v[,1]
       v.stab=FALSE
       u.stab=FALSE
       iter=0

       #--computing nx(degree of sparsity)--#
       nx = p-keepX[h]
       vect.keepX[h]=keepX[h]

       #--iterations on v and u--#
       while((v.stab==FALSE) || (u.stab==FALSE)){
            iter=iter+1
            u.old=u.new
            v.temp=t(X.temp)%*%u.old
            v.old=v.new
            
            if(h>=2){
               u.new=(lsfit(y=X%*%v.old, x=X%*%mat.v[,1:(h-1)],intercept=FALSE)$res)
               u.new=u.new/sqrt(drop(crossprod(u.new)))
            }
            
            #--penalisation on loading vectors--#
            if(nx!=0){
               v.new = ifelse(abs(v.temp) > abs(v.temp[order(abs(v.temp))][nx]), 
               (abs(v.temp) - abs(v.temp[order(abs(v.temp))][nx])) * sign(v.temp), 0)
            }
            
            if(h==1){  
               u.new = as.vector(X.temp %*% v.new)
               u.new=u.new/sqrt(drop(crossprod(u.new)))
            }
            
            #--checking convergence--#
            if(crossprod(u.new-u.old)<tol){u.stab=TRUE}
            if(crossprod(v.new-v.old)<tol){v.stab=TRUE}
            
            if ((is.na(v.stab)) | (is.na(u.stab)) | (iter >= max.iter))
            {v.stab=TRUE; u.stab=TRUE}
        
       }##fin while v


       v.final = v.new/sqrt(drop(crossprod(v.new)))
             
       #--deflation of data--#
       X.temp= X.temp - svd.X$d[1] * svd.X$u[,1] %*% t(svd.X$v[,1])
       
       
       vect.iter[h]=iter
       mat.v[,h]=v.final
       mat.u[,h]=u.new
       
       #--calculating adjusted variances explained--#
       X.var = X %*% mat.v[,1:h]%*%solve(t(mat.v[,1:h])%*%mat.v[,1:h])%*%t(mat.v[,1:h])
       vect.varX[h] = sum(X.var^2)

# KA: to add if biplot function (but to be fixed!)
       #sdev[h] = sqrt(svd.X$d[1])
       
        
    }#fin h
    
    rownames(mat.u) = ind.names

    cl = match.call()
		cl[[1]] = as.name('spca')

    result = (list(call = cl, X = X,
		   ncomp = ncomp,	
                   #sdev = sdev,  # KA: to add if biplot function (but to be fixed!)
                   #center = center, # KA: to add if biplot function (but to be fixed!)
                   #scale = scale,   # KA: to add if biplot function (but to be fixed!)
                   varX = vect.varX/sum(X^2),
                   keepX = vect.keepX,
                   iter = vect.iter,
                   rotation = mat.v,
                   x = mat.u,
                   names = list(X = X.names, indiv = ind.names),
                   loadings=list(mat.v),
                   variates=list(mat.u)
              ))
			  
    class(result) = c("spca", "prcomp", "pca")
    return(invisible(result))
}
