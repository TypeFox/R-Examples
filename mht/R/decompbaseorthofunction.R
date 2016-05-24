# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: pre 01-01-2013
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

decompbaseortho=function(data)
{
    #   data = Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.


    #		-------------------------------------
    #			checking entries
    #		-------------------------------------
    
    if(!is.matrix(data)) stop(" `data' has to be a matrix")
    if(ncol(data)<2){stop("only one dimension, no point running this algorithm.. ")}
           
    U=svd(data[,1])$u
    p=ncol(data)
    n=nrow(data)
    nonind=numeric(0)   # record the variables that don't contribute to the basis ("useless variables")
    trueind=numeric(0)  # record the variables that contribute to the basis ("useful variables")
    rank=rankMatrix(data)
    for(k in 1:(p-1))
    {
        m=0
        u=U         # orthogonal basis of Vk
        ncol_u=dim(u)[2]

        # decomposition of data(k+1) in the basis U:
        x=matrix(0,ncol_u,1)
            for(i in 1:ncol_u)
            {
                x[i,1]=sum(data[,k+1]*u[,i]) #en colonne on a la variable qui est decompose, et en ligne sa decomp dans la base
            }
        X=u%*%x # decomposition of the variable data(k+1) in the orthogonal basis of Vk (U)


        X_=data[,k+1]-X # calculation of the orthogonal part of data(k+1)
        
        if(max(abs(X_))<10^-10)		# threshold due to the calculation precision in R. We put the column to zero if under the threshold
        {
            X_=0
            u1=cbind(rep(0,nrow(data)))
            colnames(u1)="nonind"
        }else{u1=svd(X_)$u}
        
        # add the orthogonal part of data(k+1) to the orthogonal basis
        U=cbind(U,u1)
        if(sum(t(u1)%*%u1)==0){nonind=c(nonind,k+1)}else{trueind=c(trueind,k+1)}
        
        # if we reach the rank of the matrix data, no need to keep on as every other variable will be "useless"
        if((dim(U)[2]-length(nonind))==n){break}
        if(length(trueind)==(rank-1)){break}
    }
    
    # complete the matrix U with as zero columns as needed
    a=dim(U)[2]
    U=cbind(U,matrix(0,n,p-a))
    if(a<p)
    {
        nonind=c(nonind,(a+1):p)
    }
    return(list(U=U,nonind=nonind,trueind=trueind,rank=rank))
}