# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research, Toulouse France and 
# The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD 
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


pls <-
function(X, 
         Y, 
         ncomp = 2, 
         mode = c("regression", "canonical", "invariant", "classic"),
         max.iter = 500, 
         tol = 1e-06,
         near.zero.var = TRUE)
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    {
        stop("'X' must be a numeric matrix.")
    }
    
    X = as.matrix(X)
    Y = as.matrix(Y)
     
    if (!is.numeric(X) || !is.numeric(Y)) 
    {
        stop("'X' and/or 'Y' must be a numeric matrix.")
    }
    
    n = nrow(X)
    q = ncol(Y)
     
    if ((n != nrow(Y))) 
    {
        stop("unequal number of rows in 'X' and 'Y'.")
    }
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    {
        stop("invalid number of variates, 'ncomp'.")
    }
    
    if(near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position,drop=FALSE]

            if(ncol(X)==0)
            {
                stop("No more predictors after Near Zero Var has been applied!")
            }
            
        }
    }
    
	p = ncol(X)
	
    ncomp = round(ncomp)
    if(ncomp > p)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", p, ".")
        ncomp = p
    }
     
    mode = match.arg(mode)
     
     #-- initialisation des matrices --#
     X.names = dimnames(X)[[2]]
     if (is.null(X.names))
     {
         X.names = paste("X", 1:p, sep = "")
         dimnames(X)[[2]]=X.names
         
     }
     
     if (dim(Y)[2] == 1)
     {
         Y.names = "Y"
     }else{
         Y.names = dimnames(Y)[[2]]
         if (is.null(Y.names))
         {
             Y.names = paste("Y", 1:q, sep = "")
             dimnames(Y)[[2]]=Y.names
         }
     }
     
     ind.names = dimnames(X)[[1]]
     if (is.null(ind.names))
     {
         ind.names = dimnames(Y)[[1]]
         rownames(X) = ind.names
     }
     
     if (is.null(ind.names))
     {
         ind.names = 1:n
         rownames(X) = rownames(Y) = ind.names
     }		
    	
    #-- centrer et r?duire les donn?es --#
    X = scale(X, center = TRUE, scale = TRUE)
    Y = scale(Y, center = TRUE, scale = TRUE) 
     
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

      X.aux = X.temp
      if (na.X)
      {
        X.aux[is.na.X] = 0
      }

      Y.aux = Y.temp
      if (na.Y)
      {
        Y.aux[is.na.Y] = 0
      }  
      
      #-- initialisation --#
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
        # update 5.0-2: t is not normed
        #t = t / drop(sqrt(crossprod(t)))
      }else{
        t = X.temp %*% a.old / drop(crossprod(a.old))
      }
      
      if (na.Y)
      {
        u = Y.aux %*% b.old
        B = drop(b.old) %o% n.ones
        B[t(is.na.Y)] = 0
        b.norm = crossprod(B)
        u = u / diag(b.norm)
        # update 5.0-2: u is not normed
        #u = u / drop(sqrt(crossprod(u)))
      }else{
        u = Y.temp %*% b.old / drop(crossprod(b.old))
      }
      iterh = 1
         
      #-- convergence of a  --#
      repeat{
        if (na.X)
        {
          a = t(X.aux) %*% u
        }else{
          a = t(X.temp) %*% u #/ drop(crossprod(u)), useless because a is scaled after soft_thresholding
        }
        a = a / drop(sqrt(crossprod(a)))
        
        if (na.X)
        {
          t = X.aux %*% a
          A = drop(a) %o% n.ones
          A[t(is.na.X)] = 0
          a.norm = crossprod(A)
          t = t / diag(a.norm)
          # update 5.0-2: t is not normed
          #t = t / drop(sqrt(crossprod(t)))
        }else{
          t = X.temp %*% a / drop(crossprod(a))
        }
        
        if (na.Y)
        {
          b = t(Y.aux) %*% t
        }else{
          b = t(Y.temp) %*% t #/ drop(crossprod(t)), useless because b is scaled after soft_thresholding
        }
        b = b / drop(sqrt(crossprod(b)))
        
        if (na.Y)
        {
          u = Y.aux %*% b
          B = drop(b) %o% n.ones
          B[t(is.na.Y)] = 0
          b.norm = crossprod(B)
          u = u / diag(b.norm)
          # update 5.0-2: u is not normed
          #u = u / drop(sqrt(crossprod(u)))
        }else{
          u = Y.temp %*% b / drop(crossprod(b))
        }
        
        if (crossprod(a - a.old) < tol) {break}
        if (iterh == max.iter)
        {
          warning(paste("Maximum number of iterations reached for the component", h),call. = FALSE)
          break
        }
        
        a.old = a
        b.old = b
        iterh = iterh + 1
      }
        
        #-- deflation des matrices --#
        if (na.X)
        {
            X.aux = X.temp
            X.aux[is.na.X] = 0
            c = crossprod(X.aux, t)		
            T = drop(t) %o% p.ones
            T[is.na.X] = 0
            t.norm = crossprod(T)				
            c = c / diag(t.norm)
        }else{
            c = crossprod(X.temp, t) / drop(crossprod(t))
        }	
		
        X.temp = X.temp - t %*% t(c)   
         
        #-- mode canonique --#
        if (mode == "canonical")
        {
            if (na.Y)
            {
                Y.aux = Y.temp
                Y.aux[is.na.Y] = 0
                e = crossprod(Y.aux, u)
                U = drop(u) %o% q.ones
                U[is.na.Y] = 0
                u.norm = crossprod(U)				
                e = e / diag(u.norm)					
            }else{
                e = crossprod(Y.temp, u) / drop(crossprod(u))
            }
			
            Y.temp = Y.temp - u %*% t(e)
        }
         
        #-- mode classic --#
        if(mode == "classic")
        {
            Y.temp = Y.temp - t %*% t(b)
        }
        #-- mode regression --#
        if(mode == "regression")
        {
            if (na.Y)
            {
                Y.aux = Y.temp
                Y.aux[is.na.Y] = 0
                d = crossprod(Y.aux, t)
                T = drop(t) %o% q.ones
                T[is.na.Y] = 0
                t.norm = crossprod(T)				
                d = d / diag(t.norm)
            }else{
                d = crossprod(Y.temp, t) / drop(crossprod(t))
            }
             
            Y.temp = Y.temp - t %*% t(d)
        }
		
        #-- mode invariant --#
        if (mode == "invariant") {Y.temp = Y}
        
        mat.t[, h] = t
        mat.u[, h] = u
        mat.a[, h] = a
        mat.b[, h] = b
        mat.c[, h] = c
        if (mode == "regression") {mat.d[, h] = d}
	
        if (mode == "canonical") {mat.e[, h] = e}
        
        iter=c(iter,iterh) #save the number of iteration per component
    } #-- fin boucle sur h --#

    #-- valeurs sortantes --#
    rownames(mat.a) = rownames(mat.c) = X.names
    rownames(mat.b) = Y.names
    rownames(mat.t) = rownames(mat.u) = ind.names

    comp = paste("comp", 1:ncomp)
    colnames(mat.t) = colnames(mat.u) = comp
    colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp 

    cl = match.call()
    cl[[1]] = as.name('pls')
     
    result = list(call = cl,
	              X = X, 
	              Y = Y, 
	              ncomp = ncomp, 
	              mode = mode, 
	              mat.c = mat.c,
                  mat.d = mat.d,
                  mat.e = mat.e,
	              variates = list(X = mat.t, Y = mat.u),
	              loadings = list(X = mat.a, Y = mat.b), 
                  names = list(X = X.names, Y = Y.names, indiv = ind.names),
                  tol = tol,
                  max.iter = max.iter,
                  iter=iter
                )
    if (near.zero.var == TRUE) result$nzv = nzv
	
    class(result) = "pls"
    return(invisible(result))
}

