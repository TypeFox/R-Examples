# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 28-05-2014
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


# match keepX.constraint and X after removing variables
# random.subsampling(Y)
# prediction.formula(data,H,uloadings,vloadings,CH,means.X,means.Y,sigma.X,sigma.Y)
# pre.screening(X,coeff)


#match.keepX.constraint
match.keepX.constraint=function(X,names.remove,keepX.constraint)
{
    #matching the new X (after removing some variables) and keepX.constraint
    if(length(names.remove)>1)
    {
        
        ind.match= lapply(keepX.constraint,function(x){match(x,names.remove)})
        
        if(sum(!is.na(unlist(ind.match)))>0)
        {
            warnings("at least one variable was removed from keepX.constraint because of a null variance. Please check object$keepX.constraint to see which variables are used.")
            #remove from keepX.constraint
            keepX.constraint=lapply(keepX.constraint,function(x){temp=match(x,names.remove); if(sum(!is.na(temp))>0){x=x[-which(!is.na(temp))]}else{x=x}})
        }
        
        keepX=unlist(lapply(keepX.constraint,length))
        if(any(keepX==0))
        {
            ind.min=min(which(keepX==0))
            ncomp=1:(ind.min-1)
            warnings(paste("Only", ncomp,"components are used."))
            #construction of the new keepX.constraint, using ncomp components
            keepX.constraint.temp=keepX.constraint
            for(i in (ncomp+1):length(keepX.constraint)) keepX.constraint.temp[[i]]=NULL
            keepX.constraint=keepX.constraint.temp
            
        }
        
    }
    
    out=keepX.constraint
}


# random subsamplings
random.subsampling=function(Y,ratio)
{
    if(missing(ratio)) ratio=0.3

    N=length(Y)

    removed=floor(N*ratio)
    N.pick=N-removed

    f=N-signif(N.pick,2) #number of samples to discard
    keep=N-f # number of samples to keep


    a=1-round(f/N,2) # percentage of samples to keep
    n=N
    n1=sum(Y==levels(Y)[1])
    n2=sum(Y==levels(Y)[2])
    pick1=floor(a*n1)
    pick2=floor(a*n2)

    if(pick1+pick2!=keep)
    {
        a=keep-(pick1+pick2) #should be 1? or max 2? (depend on N and round)
        if(pick1>pick2) pick1=pick1+a
        if(pick2>pick1) pick2=pick2+a
    }
    #at this stage we have pick1 and pick2, number of samples to keep in each class

    a=sample(which(Y==levels(Y)[1]),replace=FALSE) #stratified bootstrap in Other
    b=sample(which(Y==levels(Y)[2]),replace=FALSE) #stratified bootstrap in MSC
    A1=a[1:pick1]
    A2=b[1:pick2]
    A=sort(c(A1,A2))

    #A contains the sample to keep for the learning set
    out=A
}



# prediction formula
prediction.formula=function(X.test,ncomp,Y.scaled,unmap.Y,variates.X,uloadings,CH,means.X,means.Y,sigma.X,sigma.Y,method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
{
    if(missing(X.test)|missing(ncomp)|missing(Y.scaled)|missing(uloadings)|missing(variates.X)|missing(CH)|missing(means.X)|missing(means.Y)|missing(sigma.X)|missing(sigma.Y)) stop("at least one argument is missing")

    if(missing(method)) method="all"
    betay = list()
    nlevelY=ncol(unmap.Y)
    B.hat=array(0,c(ncol(X.test),nlevelY,ncomp))
    Y.hat=array(0,c(nrow(X.test),nlevelY,ncomp))
    q=nlevelY #number of columns of Y, just 2 here
    t.hat=array(0,c(nrow(X.test),ncomp))

    for(num.comp in 1:ncomp){
        
        dd= coefficients(lm(Y.scaled~variates.X[,1:num.comp,drop=FALSE])) #regression of Y on variates.global.X => =loadings.global.Y at a scale factor
        if(q>=2){betay[[num.comp]]=(dd[-1,])}
        
        W = uloadings[, 1:num.comp,drop=FALSE] %*% solve(t(CH[, 1:num.comp,drop=FALSE]) %*% uloadings[, 1:num.comp,drop=FALSE])
        B = W %*% drop(betay[[num.comp]])
        
        Y.temp=scale(X.test,center=means.X,scale=sigma.X) %*% as.matrix(B) #so far: gives a prediction of Y centered and scaled
        Y.temp2=scale(Y.temp,center=FALSE,scale=1/sigma.Y) #so far: gives a prediction of Y centered, with the right scaling
        Y.temp3=scale(Y.temp2,center=-means.Y,scale=FALSE) #so far: gives a prediction of Y with the right centering and scaling
    
        
        Y.hat[, , num.comp] = Y.temp3 # we add the variance and the mean of Y used in object to predict
        t.hat[, num.comp] = scale(X.test, center = means.X, scale = sigma.X) %*% W[, num.comp]
        B.hat[, , num.comp] = B
    }  #end h

#predicted=apply(Y.hat,c(1,3),which.max)

    G = matrix(0, nrow = q, ncol = ncomp)
    cls = list()

    for (i in 1:q) {
        if(ncomp > 1) {
            G[i, ] = apply(variates.X[unmap.Y[, i] == 1, ], 2, mean)
        }
        else {
            G[i, ] = mean(variates.X[unmap.Y[, i] == 1, ])
        }
    }

    # ----    max distance -----------------

    if (any(method == "all") || any(method == "max.dist")) {
        
        function.pred = function(x){
            nr = nrow(x)
            tmp = vector("numeric", nr)
            for(j in 1:nr){
                tmp[j] = (which(x[j, ] == max(x[j, ]))[1])
            }
            return(tmp)
        }
        cls$max.dist = matrix(apply(Y.hat, 3, function.pred), ncol = ncomp)
        colnames(cls$max.dist) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
    }

    # ----    centroids distance -----------------

    if (any(method == "all") || any(method == "centroids.dist")) {
        
        cl = matrix(nrow = nrow(X.test), ncol = ncomp)
        
        centroids.fun = function(x, G, h) {
            q = nrow(G)
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                d = apply((x - G[, 1:h])^2, 1, sum)
            }
            else {
                d = (x - G[, 1])^2
            }
            cl.id = which.min(d)
        }
        
        for (h in 1:ncomp) {
            cl.id = apply(matrix(t.hat[, 1:h], ncol = h), 1, centroids.fun, G = G, h = h)
            cl[, h] = cl.id
        }
        colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
        cls$centroids.dist = cl
    }

    # ----    mahalanobis distance -----------------

    if (any(method == "all") || any(method == "mahalanobis.dist")) {
        
        cl = matrix(nrow = nrow(X.test), ncol = ncomp)
        
        Sr.fun = function(x, G, unmap.Y, h) {
            q = nrow(G)
            Xe = unmap.Y %*% G[, 1:h]
            Xr = variates.X[, 1:h] - Xe
            Sr = t(Xr) %*% Xr / nrow(unmap.Y)
            Sr.inv = solve(Sr)
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                mat = (x - G[, 1:h]) %*% Sr.inv %*% t(x - G[, 1:h])
                d = apply(mat^2, 1, sum)
            }
            else {
                d = drop(Sr.inv) * (x - G[, 1])^2
            }
            cl.id = which.min(d)
        }
        
        for (h in 1:ncomp) {
            cl.id = apply(matrix(t.hat[, 1:h], ncol = h), 1, Sr.fun, G = G, unmap.Y = unmap.Y, h = h)
            cl[, h] = cl.id
        }
        colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
        cls$mahalanobis.dist = cl
    }
    #-- valeurs sortantes --#
    if (any(method == "all")) method = "all"
    rownames(t.hat) = rownames(X.test)
    colnames(t.hat) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(X.test)
    colnames(Y.hat) = colnames(unmap.Y)
    dimnames(Y.hat)[[3]]=paste("comp.",1:ncomp,sep="")
    colnames(G) = paste("dim", c(1:ncomp), sep = " ")
    
    out=list(Y.hat=Y.hat,t.hat=t.hat,B.hat=B.hat,centroids = G,method = method, class = cls)

}
