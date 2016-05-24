# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD 
# created: pre 04-08-2014
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



# -------------------------------------------------------------------------------
# helpers for pls and spls
# -------------------------------------------------------------------------------

# --------------------------------------
# Check.entry.bootsPLS
# --------------------------------------

Check.entry.bootsPLS = function(X, Y)
{


    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")

    X = as.matrix(X)

    if (!is.numeric(X))
    stop("'X' must be a numeric matrix.")

    if(!is.factor(Y)) stop("Y must be a factor")


    N = nrow(X)
    P= ncol(X)

    if ((N != length(Y)))
    stop("unequal number of rows in 'X' and 'Y'.")


    #-- initialisation des matrices --#
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]]=X.names
    }


    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = dimnames(Y)[[1]]
        rownames(X) = ind.names
    }

    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X) = rownames(Y) = ind.names
    }

    if(length(unique(X.names))!=P) stop("Unique indentifier is needed for the columns of X")

    return(list(X=X,Y=Y,X.names=X.names,ind.names=ind.names))


}

# --------------------------------------
# Check.entry.X
# --------------------------------------
Check.entry.X = function(X)
{
    
    
    #  if(length(levels(study)) == 1)  # Aida
    #  stop("\nstudys must have more than one level")      #WHY?
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
   
    N = nrow(X)
    P= ncol(X)


    #-- initialisation des matrices --#
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]]=X.names
    }
    
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X) = ind.names
    }
    
    if(length(unique(X.names))!=P) stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X,X.names=X.names,ind.names=ind.names))
}


# --------------------------------------
# Check.entry.pls
# --------------------------------------

Check.entry.pls = function(X, Y, ncomp, keepX, keepY)
{
    
    
    #  if(length(levels(study)) == 1)  # Aida
    #  stop("\nstudys must have more than one level")      #WHY?
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    Y = as.matrix(Y)
    
    if (!is.numeric(X) || !is.numeric(Y))
    stop("'X' and/or 'Y' must be a numeric matrix.")
    
    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)

    if ((N != nrow(Y)))
    stop("unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    if (length(keepX) != ncomp)
    stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    
    if (length(keepY) != ncomp)
    stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepY > ncol(Y)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    

     ncomp = round(ncomp)
     if(ncomp > P)
     {
         warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
         ncomp = P
     }
    
    
     #-- initialisation des matrices --#
     X.names = dimnames(X)[[2]]
     if (is.null(X.names))
     {
         X.names = paste("X", 1:P, sep = "")
         dimnames(X)[[2]]=X.names
     }

     if (dim(Y)[2] == 1) Y.names = "Y"
     if (dim(Y)[2] > 1)
     {
         Y.names = dimnames(Y)[[2]]
         if (is.null(Y.names))
         {
             Y.names = paste("Y", 1:Q, sep = "")
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
         ind.names = 1:N
         rownames(X) = rownames(Y) = ind.names
     }
     
     if(length(unique(X.names))!=P) stop("Unique indentifier is needed for the columns of X")
     if(length(unique(Y.names))!=Q) stop("Unique indentifier is needed for the columns of Y")

    return(list(X=X,Y=Y,ncomp=ncomp,X.names=X.names,Y.names=Y.names,ind.names=ind.names))
}



# --------------------------------------
# l2.norm
# --------------------------------------
l2.norm=function(x)
{
    if(!is.vector(x)) stop("x has to be a vector")
    out=x/drop(sqrt(crossprod(x)))
}



# --------------------------------------
# soft_thresholding
# --------------------------------------
soft_thresholding=function(x,nx)
{
    #if (nx != 0) {
    #    x = ifelse(abs(x) > abs(x[order(abs(x))][nx]),
    #    (abs(x) - abs(x[order(abs(x))][nx])) * sign(x), 0)
    #}
    
    #selection on a (loadings.X). modified on 19/02/15 to make sure that a!=0
    if(nx!=0)
    {
        absx=abs(x)
        if(sum(rank(absx)<=nx)>0)# if nx is not high enough, we don't put any coefficients to zero
        {
            x=ifelse(absx>absx[which(rank(absx)==max(rank(absx)[which(rank(absx)<=(nx))]))[1]],
            (absx-absx[which(rank(absx)==max(rank(absx)[which(rank(absx)<=(nx))]))[1]])*sign(x),0)
        }
    }
    
    x
}

