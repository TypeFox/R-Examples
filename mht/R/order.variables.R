# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: pre 01-01-2013
# last modification: 14-03-2015
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

order.variables=function(data,Y,maxordre,ordre=c("bolasso","pval","pval_hd","FR"),var_nonselect,m,showordre)
{
	
    #-----------------------------------	
    #   data = Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.
    #   Y = Response variable of length n.
    #	maxordre= Number of variables to be ordered in the first step of the procedure
    #	ordre= which method should be used to order the variables. One of pval, pval_hd, FR or bolasso
    #	var_nonselect= Number of variables that are not submitted to feature selection; they are the first columns of the data matrix
    #	m= number of replications of the Lasso in the Bolasso technique. Only used when ordre is set to `bolasso'
    #	showordre=affiche l'ordre au fur et a mesure
    #-----------------------------------------

    data = as.matrix(data)
    n=nrow(data)
    p=ncol(data)
    
    Y=as.numeric(Y)
    #-- validation des arguments --#
    if (length(dim(data)) != 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix.")
    if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}

    #		-------------------------------------
    #			rank the variables
    #		-------------------------------------

    if(missing(m)){m=100}
    if(missing(maxordre)){maxordre=min(n/2-1,p/2-1)}
    if(missing(ordre)){ordre="bolasso"}
    if(missing(showordre)){showordre=1}

    

    #		-------------------------------------
    #			scaling the data matrix
    #		-------------------------------------
    
    #si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute et on rajoute 1 dans var_nonselect s'il n'etait pas manquant
    temp=data.scale(data)
    data=temp$data
    intercept=temp$intercept
    means.X=temp$means.data
    sigma.X=temp$sigma.data
    p=ncol(data)
    
    if(missing(var_nonselect))
    {
        var_nonselect=1
    }else{
        if(!intercept){var_nonselect=var_nonselect+1}
    }
    if(var_nonselect<0){stop("var_nonselect has to be nonnegative")}
    
    maxordre=min(n-1,p-1,maxordre+!intercept)


    if(ordre=="pval")
    {
        if(p<n) 	#on calcule pval avec une seule regression comportant toutes les variables
        {reg=lm(Y~data-1)
            a=summary(reg)$coefficients[,4]
            a[1:var_nonselect]=0	# no feature selection on the var_nonselect first variables
            b=sort(a,index.return=TRUE)
            ORDRE=b$ix[1:maxordre]
            ORDREBETA=matrix(b$ix,nrow=1)
            data_ord=data[,b$ix] #on a ainsi les XI ordonnees
            if(showordre){print(colnames(data)[b$ix])}
            }else{		#on calcule pval avec FDR2: une regression pour chaque variable

                print("p>n, the order 'pval' is not possible, 'pval_hd' is used instead")
                ordre="pval_hd"
            }
    }
    if(ordre=="pval_hd")
    {
        a=numeric(0)
        for(i in 1:p)
        {
            reg=lm(Y~data[,i]-1)
            c=summary(reg)$coefficients[,4]
            a=c(a,c)
        }
        a[1:var_nonselect]=0	# no feature selection on the var_nonselect first variables
        b=sort(a,index.return=TRUE)
        ORDREBETA=matrix(b$ix,nrow=1)
        ORDRE=b$ix[1:maxordre]
        data_ord=data[,b$ix] #on a ainsi les XI ordonnees
        if(showordre){print(colnames(data)[b$ix])}
    }


    if(ordre=="bolasso")
    {
        res=dyadiqueordre(data,Y,m,maxordre,var_nonselect,showordre=showordre)# donne l'ordre (dans ordre) et le nombre de fois ou l'algo a redemarre (dans prob)
        b=res$ordre
        ORDRE=b

        a=match(1:p,b)
        b=c(b,(1:p)[which(is.na(a))])

        ORDREBETA=matrix(b,nrow=1)
        data_ord=data[,b] #on a ainsi les XI ordonnees dans data_ord
    }

    if(ordre=="FR")
    {
        path=1
        ind.left=2:p
        while(length(path)<min(n-1,p,maxordre))
        {
            res=NULL
            for(i in 1:length(ind.left))
            {
                X=data[,c(path,ind.left[i])]
                I=diag(1,n)
                H=X%*%solve(t(X)%*%X)%*%t(X)
                res[i]=t(Y)%*%(I-H)%*%Y	#same as reg=lm(Y~XI[,VRAINDICE]-1); sum(reg$residuals^2)
            }
            a=which.min(res)
            path=c(path,ind.left[a])
            ind.left=ind.left[-a]	
            
        }
        b=ORDRE=path
        a=match(1:p,b)
        b=c(b,(1:p)[which(is.na(a))])
        
        ORDREBETA=matrix(b,nrow=1)
        data_ord=data[,b] #on a ainsi les XI ordonnees dans data_ord
        if(showordre){print(colnames(data)[path])}
    }

#ORDRE is the order on the maxordre first variables of data
#ORDREBETA is the order on all the variables of data (either arbitrary completion of ORDRE, or the true order with pval, pval_hd)
#data_ord is the data ordered by ORDREBETA

return(list(data=list(X=data,Y=Y,means.X=means.X,sigma.X=sigma.X),data_ord=data_ord,ORDRE=ORDRE,ORDREBETA=ORDREBETA))
}
