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

bolasso<-function(data,Y,mu,m,probaseuil,penalty.factor,random)
{
    #-----------------------------------
    #   data = Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.
    #   Y = Response variable of length n.
    #   mu = Positive regularization sequence to be used for the Lasso.
    #   m = Number of bootstrapped iteration of the Lasso. Default is m=100.
    #   probaseuil = A frequency threshold for selecting the set of relevant variables. Default is 1.
    #   penalty.factor = Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables except the intercept.
    #   random = optionnal parameter, matrix of size n*m. If \code{random} is provided, the m bootstrap samples are constructed from its m columns.
    #-----------------------------------------


    #		-------------------------------------
    #			checking entries
    #		-------------------------------------
    
    if(!is.matrix(data)) stop(" `data' has to be a matrix")
    if(is.matrix(Y)){
        if(ncol(Y)>1) stop("`Y' has to be a vector or a single column matrix")
    }
    ntot=nrow(data)
    p=ncol(data)
    if(length(Y)!=ntot){stop(" 'data' and 'Y' must have the same length ")}

    if(missing(mu)) stop("`mu' is missing")
    if(missing(probaseuil)){probaseuil=1}
    if(missing(m)){m=100}
    
    if(any(mu<0)) stop(" `mu' has to be a positive regularization sequence")
    if(m<1) stop(" `m' has to be greater than 1")
    if( probaseuil<0 | probaseuil>1) stop("`probaseuil' has to be in (0,1)")
    if(!missing(random))
    {
        if(!is.matrix(random)) stop("`random' has to be a matrix")
        if(ncol(random)!=m) stop("`random' has to have `m' columns")
    }
   
    


    #		-------------------------------------
    #			scaling `data', adding intercept if needed
    #		-------------------------------------
    temp=data.scale(data)
    data=temp$data
    intercept=temp$intercept
    means.X=temp$means.data
    sigma.X=temp$sigma.data
    p=ncol(data)
                
    if(missing(penalty.factor))
    {
        penalty.factor=c(0,rep(1,p-1))
    }else{
        if(length(penalty.factor)>p) stop("`penalty.factor' has to be a vector of length ncol(data)")
        if(!intercept){penalty.factor=c(0,penalty.factor)}
    }
        

    dimu=length(mu)
    # first fit of the model with Lasso
    lasso1=glmnet(data,Y,alpha=1,lambda=mu,penalty.factor=penalty.factor)
    mu=lasso1$lambda # gives the same order for the regularization parameter between glmnet and these outputs

    mat=NULL

    # for each mu -penalization/regularization parameter, we fit a linear model on the selected coefficients (from the lasso at penalty `mu') and we record the fitted values in `mat'
    for(i in 1:dimu)
    {
        if(penalty.factor[1]==0)
        {
            beta_ind=c(1,which(lasso1$beta[,i]!=0))
        }else{
            beta_ind=which(lasso1$beta[,i]!=0)
        }

        reg=lm(Y~data[,beta_ind]-1)
        beta0=reg$coefficients
        beta0[-which(beta0!=0)]=0
        mat0=as.matrix(data[,beta_ind])%*%beta0
        mat=cbind(mat,mat0)
    }
    
    # we record the residuals of the model for each `mu'
    eps0=as.matrix(Y)%*%matrix(1,1,dimu)-mat
    eps1=scale(eps0,center=TRUE,scale=FALSE)

    # if the matrix `random' was not supplied, we construct a random matrix from which we will perform the subsamplings
    if(missing(random))
    {
        a=matrix(runif(ntot*m,0,ntot),nrow=ntot)
        random=ceiling(a)
    }


    # for the `m' subsamplings, the lasso is performed on `data 'and a bootstrapped Y (on residuals), the selected variables are recorded in `compteur'
    compteur=matrix(0,p,dimu)
    for (j in 1:m)
    {
        Y_boot=as.matrix(mat+eps1[random[,j],])	#on bootstrap les residus
        
        betaboot2=matrix(0,p,dimu)
        for(i in 1:dimu)
        {
            mu1=mu[i]
            lasso1=glmnet(data,Y_boot[,i],alpha=1,lambda=mu1,penalty.factor=penalty.factor)
        
            betaboot2[,i]=(lasso1$beta[,1]!=0)
            betaboot2[which(penalty.factor==0),i]=1
        }
        compteur=compteur+betaboot2 #record of the selected variables

    }


    compteur2=as.matrix(compteur)

    probavariable=compteur2/m
    colnames(probavariable)=mu
    beta_ind=(probavariable>=probaseuil)
    rownames(probavariable)=rownames(beta_ind)=colnames(data)

    out=list(data=list(X=data,Y=Y,means.X=means.X,sigma.X=sigma.X),ind=beta_ind,frequency=probavariable,mu=mu,call=match.call())
    out
    structure(out,class="bolasso")
    
}