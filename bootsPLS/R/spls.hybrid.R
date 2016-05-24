# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 28-05-2014
# last modification: 16-03-2015
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


# perform a ncomp-1 pls on the keepX.constraint variables to keep, and a spls with keepX on the last component
spls.hybrid <-function(X,
                Y,
                ncomp = 2,
                mode = c("regression", "canonical", "invariant", "classic"),
                max.iter = 500,
                tol = 1e-06,
                keepX.constraint,
                keepY.constraint, 
                keepX,
                keepY,
                near.zero.var = FALSE)
{
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    Y = as.matrix(Y)
    
    if (!is.numeric(X) || !is.numeric(Y))
    stop("'X' and/or 'Y' must be a numeric matrix.")
    
    
    if(missing(keepX.constraint))
    {
        if(missing(keepX))
        {
            keepX=rep(ncol(X),ncomp)
        }
        keepX.constraint=list()
    }else{
        if(missing(keepX))
        {
            keepX=NULL
        }
    }
    
    if(missing(keepY.constraint))
    {
        if(missing(keepY))
        {
            keepY=rep(ncol(Y),ncomp)
        }
        keepY.constraint=list()
    }else{
        if(missing(keepY))
        {
            keepY=NULL
        }
    }
    
    
    if((length(keepX.constraint)+length(keepX))!=ncomp) stop("length (keepX.constraint) + length(keepX) should be ncomp")
    if((length(keepY.constraint)+length(keepY))!=ncomp) stop("length (keepY.constraint) + length(keepY) should be ncomp")
    
    keepX.temp=c(unlist(lapply(keepX.constraint,length)),keepX) #of length ncomp
    keepY.temp=c(unlist(lapply(keepY.constraint,length)),keepY) #of length ncomp
    
    
    #-- validation des arguments --#
    check=Check.entry.pls(X,Y,ncomp,keepX.temp,keepY.temp) # gives colnames and rownames to X, Y if none at first
    X=check$X
    Y=check$Y
    ncomp=check$ncomp
    X.names=check$X.names
    Y.names=check$Y.names
    ind.names=check$ind.names
    
    # match keepX.constraint and the colnames of X in order for keepX.constraint to be a list of character
    if(length(keepX.constraint)>0)
    {
        X.indice=X[,unlist(keepX.constraint),drop=FALSE]
        keepX.constraint=relist(colnames(X.indice),skeleton=keepX.constraint)
    }
    
    # same for keepY.constraint
    if(length(keepY.constraint)>0)
    {
        Y.indice=Y[,unlist(keepY.constraint),drop=FALSE]
        keepY.constraint=relist(colnames(Y.indice),skeleton=keepY.constraint)
    }
    
    if(near.zero.var == TRUE)
    {
        nzv.X = nearZeroVar(X)
        if (length(nzv.X$Position > 0))
        {
            names.remove.X=colnames(X)[nzv.X$Position]
            warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv$X for problematic predictors.")
            X = X[, -nzv.X$Position]
            if(ncol(X)==0) {stop("No more variables in X")}
        }else{names.remove.X=NULL}
        
        nzv.Y = nearZeroVar(Y)
        if (length(nzv.Y$Position > 0))
        {
            names.remove.Y=colnames(Y)[nzv.Y$Position]
            warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv$Y for problematic predictors.")
            Y = Y[, -nzv.Y$Position]
            if(ncol(Y)==0) {stop("No more variables in Y")}
        }else{names.remove.Y=NULL}
        
    }else{
        # remove columns with no variance in X. no need to remove them from keepX.constraint because it is matched to the colnames of X.
        # However, need to check that there are still variables in each keepX.constraint[[i]]
        X.scale=scale(X)
        sigma.X=attr(X.scale,"scaled:scale")
        remove=which(sigma.X==0)
        if(length(remove)>0)
        {
            names.remove.X=colnames(X)[remove]
            X=X[,-remove,drop=FALSE]
            X.names=X.names[-remove]
        }else{names.remove.X=NULL}
        if(ncol(X)==0) {stop("No more variables in X")}
        
        # remove columns with no variance in Y. no need to remove them from keepX.constraint because it is matched to the colnames of X.
        # However, need to check that there are still variables in each keepX.constraint[[i]]
        Y.scale=scale(Y)
        sigma.Y=attr(Y.scale,"scaled:scale")
        remove=which(sigma.Y==0)
        if(length(remove)>0)
        {
            names.remove.Y=colnames(Y)[remove]
            Y=Y[,-remove,drop=FALSE]
            Y.names=Y.names[-remove]
        }else{names.remove.Y=NULL}
        if(ncol(Y)==0) {stop("No more variables in Y")}
        
    }
    
    keepX.constraint=match.keepX.constraint(X,names.remove.X,keepX.constraint)
    keepY.constraint=match.keepX.constraint(Y,names.remove.Y,keepY.constraint)
    
    # we need numbers in keepX.constraint from now on
    keepX.constraint= lapply(keepX.constraint,function(x){match(x,colnames(X))})
    keepY.constraint= lapply(keepY.constraint,function(x){match(x,colnames(Y))})
    
    
    # update keepX and keepY
    keepX=c(unlist(lapply(keepX.constraint,length)),keepX) #of length ncomp, can contains 0
    keepY=c(unlist(lapply(keepY.constraint,length)),keepY) #of length ncomp, can contains 0
    
    
    n = nrow(X)
    q = ncol(Y)
    p = ncol(X)
    
    
    mode = match.arg(mode)
    
    #-- center and scale data --#
    X = scale(X, center = TRUE, scale = TRUE)
    Y = scale(Y, center = TRUE, scale = TRUE)
    means.X=attr(X,"scaled:center")
    sigma.X=attr(X,"scaled:scale")
    means.Y=attr(Y,"scaled:center")
    sigma.Y=attr(Y,"scaled:scale")
    
    X.temp = X
    Y.temp = Y
    mat.t = matrix(nrow = n, ncol = ncomp)
    mat.u = matrix(nrow = n, ncol = ncomp)
    mat.a = matrix(nrow = p, ncol = ncomp)
    mat.b = matrix(nrow = q, ncol = ncomp)
    mat.c = matrix(nrow = p, ncol = ncomp)
    mat.d = matrix(nrow = q, ncol = ncomp)
    mat.e = matrix(nrow = q, ncol = ncomp)
    n.ones = rep(1, n)
    p.ones = rep(1, p)
    q.ones = rep(1, q)
    na.X = FALSE
    na.Y = FALSE
    is.na.X = is.na(X)
    is.na.Y = is.na(Y)
    if (any(is.na.X)) na.X = TRUE
    if (any(is.na.Y)) na.Y = TRUE
    
    iter=NULL
    #-- boucle sur h --#
    for (h in 1:ncomp) {
        
        nx=p-keepX[h]
        ny=q-keepY[h]
        
        
        #-- svd de M = t(X)*Y --#
        X.aux = X.temp
        if (na.X) {X.aux[is.na.X] = 0}
        
        Y.aux = Y.temp
        if (na.Y) {Y.aux[is.na.Y] = 0}
        
        M = crossprod(X.aux, Y.aux)
        svd.M = svd(M, nu = 1, nv = 1)
        a.old = svd.M$u
        b.old = svd.M$v
        
        #-- latent variables --#
        if (na.X)
        {
            t = X.aux %*% a.old
            A = drop(a.old) %o% n.ones
            A[t(is.na.X)] = 0
            a.norm = crossprod(A)
            t = t / diag(a.norm)
        }else{
            t = X.aux %*% a.old / drop(crossprod(a.old))
        }
        
        if (na.Y)
        {
            u = Y.aux %*% b.old
            B = drop(b.old) %o% n.ones
            B[t(is.na.Y)] = 0
            b.norm = crossprod(B)
            u = u / diag(b.norm)
        }else{
            u = Y.aux %*% b.old / drop(crossprod(b.old))
        }
        
        iterh = 1
        
        #-- boucle jusqu'? convergence de a et de b --#
        repeat {
            a = t(X.aux) %*% u
            b = t(Y.aux) %*% t
            
            if(h<=length(keepX.constraint))
            {
                #constraint on a, associated to X
                if (nx != 0){a[-keepX.constraint[[h]]]=0}
                a=l2.norm(as.vector(a))
            }
            
            if(h<=length(keepY.constraint))
            {
                #constraint on b, associated to Y
                if (ny != 0){a[-keepY.constraint[[h]]]=0}
                b=l2.norm(as.vector(b))
            }
            
            if(h>length(keepX.constraint))
            {
                #penalisation on a, associated to X
                if (nx != 0){a=soft_thresholding(a,nx)}
                a=l2.norm(as.vector(a))
            }
            if(h>length(keepY.constraint))
            {
                #penalisation on b, associated to Y
                if (ny != 0){b=soft_thresholding(b,ny)}
                b=l2.norm(as.vector(b))
            }
            
            if (na.X)
            {
                t = X.aux %*% a
                A = drop(a) %o% n.ones
                A[t(is.na.X)] = 0
                a.norm = crossprod(A)
                t = t / diag(a.norm)
            }else{
                t = X.aux %*% a / drop(crossprod(a))
            }
            
            if (na.Y)
            {
                u = Y.aux %*% b
                B = drop(b) %o% n.ones
                B[t(is.na.Y)] = 0
                b.norm = crossprod(B)
                u = u / diag(b.norm)
            }else{
                u = Y.aux %*% b / drop(crossprod(b))
            }
            
            if (crossprod(a - a.old) < tol) {break}
            
            if (iterh == max.iter)
            {
                warning(paste("Maximum number of iterations reached for the component", h),
                call. = FALSE)
                break
            }
            
            a.old = a
            b.old = b
            iterh = iterh + 1
        }
        
        #-- deflation des matrices --#
        if (na.X)
        {
            c = crossprod(X.aux, t)
            T = drop(t) %o% p.ones
            T[is.na.X] = 0
            t.norm = crossprod(T)
            c = c / diag(t.norm)
        }else{
            c = crossprod(X.aux, t) / drop(crossprod(t))
        }
        
        X.temp = X.temp - t %*% t(c)
        
        #-- mode canonique --#
        if (mode == "canonical")
        {
            if (na.Y)
            {
                e = crossprod(Y.aux, u)
                U = drop(u) %o% q.ones
                U[is.na.Y] = 0
                u.norm = crossprod(U)
                e = e / diag(u.norm)
            }else{
                e = crossprod(Y.aux, u) / drop(crossprod(u))
            }
            
            Y.temp = Y.temp - u %*% t(e)
        }
        
        #-- mode regression --#
        if(mode == "regression")
        {
            if (na.Y)
            {
                d = crossprod(Y.aux, t)
                T = drop(t) %o% q.ones
                T[is.na.Y] = 0
                t.norm = crossprod(T)
                d = d / diag(t.norm)
            }else{
                d = crossprod(Y.aux, t) / drop(crossprod(t))
            }
            
            Y.temp = Y.temp - t %*% t(d)
        }
        
        mat.t[, h] = t
        mat.u[, h] = u
        mat.a[, h] = a
        mat.b[, h] = b
        mat.c[, h] = c
        if (mode == "regression") {mat.d[, h] = d}
        if (mode == "canonical") {mat.e[, h] = e}
        
        #-- mode classic --#
        if(mode == "classic") {Y.temp = Y.temp - t %*% t(b)}
        
        #-- mode invariant --#
        if (mode == "invariant") {Y.temp = Y}
        
        iter=c(iter,iterh) #save the number of iteration per component
    } #-- fin boucle sur h --#
    
    #-- valeurs sortantes --#
    rownames(mat.a) = rownames(mat.c) = X.names
    rownames(mat.b) = Y.names
    rownames(mat.t) = rownames(mat.u) = ind.names
    
    
    
    dim = paste("comp", 1:ncomp)
    colnames(mat.t) = colnames(mat.u) = dim
    colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = dim
    
    cl = match.call()
    cl[[1]] = as.name('spls')
    
    result = list(call = cl,
    X = X, Y = Y, ncomp = ncomp, mode = mode,
    keepX.constraint = keepX.constraint,
    keepY.constraint = keepY.constraint,
    keepX=keepX,
    keepY=keepY,
    mat.c = mat.c,
    mat.d = mat.d,
    mat.e = mat.e,
    variates = list(X = mat.t, Y = mat.u),
    loadings = list(X = mat.a, Y = mat.b),
    names = list(X = X.names, Y = Y.names, indiv = ind.names),
    tol = tol,
    max.iter = max.iter,iter=iter
    )
    
    
    if (near.zero.var == TRUE)
    {
        result$nzv$X = nzv.X
        result$nzv$Y = nzv.Y
    }
    result$coeff=list(means.X=means.X,sigma.X=sigma.X,means.Y=means.Y,sigma.Y=sigma.Y)
    
    class(result) = c("spls.hybrid","pls")
    return(invisible(result))
}# end function


