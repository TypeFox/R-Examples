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

mht=function(data,Y,var_nonselect,alpha,sigma,maxordre,ordre=c("bolasso","pval","pval_hd","FR"),m,show,IT,maxq)
{
	
    #-----------------------------------	
    #   data = Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.
    #   Y = Response variable of length n.
    #	var_nonselect= Number of variables that are not submitted to feature selection; they are the first columns of the data matrix
    #	alpha= Type I error of the tests
    # 	sigma= positive number if the variance is known, otherwise 0. The statistics are different for known and unknown variance
    #	maxordre= Number of variables to be ordered in the first step of the procedure
    #	ordre= which method should be used to order the variables. One of pval, pval_hd, FR or bolasso
    #	m= number of replications of the Lasso in the Bolasso technique. Only used when ordre is set to `bolasso'
    #	show=c(showordre,showresult).
        #	showordre=show the rank variables
        #	showresult=show the results of the multiple testing procedure
    #	IT=Number of simulations to calculate the quantiles
    #	maxq= number of alternative hypotheses that are to be tested
    #-----------------------------------------

    data = as.matrix(data)
    Y=as.numeric(Y)

    #-- validation des arguments --#
    if (length(dim(data)) != 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix.")

    n=nrow(data)
    p=ncol(data)
    #safety
    if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}
    if(missing(sigma)){sigma=0}else{if (sigma<0) stop("sigma is a positive quantity")}
    
    if(p<3) stop("The function cannot run for less than 3 variables")

    #set the default values
    if(missing(m)){m=100}
    if(missing(maxordre)){maxordre=min(n/2-1,p/2-1)}
    if(missing(alpha)){alpha=c(0.1,0.05)}
    if(missing(ordre)){ordre="bolasso"}
    if(missing(IT)){IT=1000}
    if(missing(maxq)){maxq=log(min(n,p)-1,2)}
    if(missing(show)){show=c(1,0,1)}
    showordre=show[1]
    showtest=show[2]
    showresult=show[3]

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
    maxq=min(maxq,log(min(n,p)-1,2)+1)



    #		-------------------------------------
    #			rank the variables
    #		-------------------------------------
    message("Step #1: ordering the variables")
    res=order.variables(data,Y,maxordre,ordre,var_nonselect,m,showordre)
    XI_ord=res$data_ord
    ORDRE=res$ORDRE
    ORDREBETA=res$ORDREBETA

    #		------------------------------------------
    #		  get an orthogonal family out of XI_ord
    #		------------------------------------------
    dec=decompbaseortho(XI_ord)# gives a orthogonal basis of  X(1),..X(p) where U[,i] is the decomposition of X(i) in the basis (Gram-Shmidt algorithm)
    nonind=dec$nonind #ex: nonind=4 means that the variable X(4) is a linear combination of X(1),X(2),X(3).
    U=dec$U     #U is the orthogonal basis of X(1),...,X(p).

    if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}
    dim_X=ncol(Uchap)			#nombre de variables utiles =dec$rank
    
    
    beta=lm(Y~Uchap-1)$coefficients # get the coefficients of the linear regression of Y onto the orthogonal basis U
    beta[-which(beta!=0)]=0

    #		----------------------------------------------------------------------------
    #			start of the sequential testing procedure
    #		----------------------------------------------------------------------------
    message("Step #2: multiple hypotheses testing")

    maxqdep=min(log(min(n,dim_X)-1,2)+1,maxq) #number of alternative hypotheses testing at the start of the procedure. Since the procedure is sequential, the number of alternative may be reduced when ktest is high
    maxq=maxqdep
    quantile=array(0,c(length(alpha),maxq,dim_X)) #record the quantile
    NBR=numeric(0) #record the number of hypotheses actually tested
    NBR_effect=numeric(0) #record the number of variables selected
    relevant_variables=numeric(0) #record the relevant variables
    indice=var_nonselect-1-sum(nonind<=var_nonselect) #on fait les calculs de quantiles pour tous les alpha en meme temps, donc pas besoin de le refaire pour chaque avec indice
    for(alph in 1:length(alpha)) #boucle sur alpha
    {
        ktest=var_nonselect-sum(nonind<=var_nonselect) # initialization with the first variables to test
        
        nbr_test=var_nonselect

        T=1
        while((T>0)&&(dim_X>ktest+1))
            {
                if(showresult)
                {
                    cat(rep("==",10),"\n")
                    cat("ktest=",ktest,"alpha=",alpha[alph],"\n")
                    cat(rep("==",10),"\n")
                }
                maxq=min(log(min(n,dim_X)-ktest-1,2)+1,maxqdep)
                
                if(ktest>indice)
                {
                    quant=quantilemht(XI_ord,k=ktest,alpha,IT,maxq=maxq,sigma=sigma) # calculation of the quantile, call a C code
                    quantile[,1:maxq,ktest+(var_nonselect==0)]=quant$quantile # store the results
                    indice=indice+1
                }
                if(quant$nbrprob==3){break}

                #test S, record the results of the maxq multiple tests in bb
                bb=numeric(0)
                statistics.quantile=NULL
                for(m in 0:(maxq-1))
                {
                    if(sigma==0)
                    {
                        a=(n-ktest-2^m)*sum(beta[(ktest+1):(ktest+2^m)]^2)/(2^m*sum((Y-as.matrix(Uchap[,1:(ktest+2^m)])%*%beta[1:(ktest+2^m)])^2))
                    }else{
                            a=sum(beta[(ktest+1):(ktest+2^m)]^2)/sigma
                    }
                    b=a>quantile[alph,m+1,ktest+(var_nonselect==0)]#F #1 if the null hypothesis is rejected, 0 otherwise
                    bb=c(bb,b) #record all the multiple tests                    }
                    statistics.quantile=cbind(statistics.quantile,c(a,quantile[alph,m+1,ktest+(var_nonselect==0)]))
                }
                statistics.quantile=rbind(statistics.quantile,bb)
                rownames(statistics.quantile)=c("Statistics","Quantile","Rejected")
                colnames(statistics.quantile)=paste("S_{(",ktest,"),(",0:(maxq-1),")}",sep="")
                

                if(showresult){print(statistics.quantile)}
            
                # if there is at least one 1, the null hypothesis is rejected, k=k+1 and we start again
                if(sum(bb)>0)
                {
                    T=1
                    ktest=ktest+1
                        
                    nbr_test=nbr_test+1
                    if(length(nonind)>0)
                    {
                        for(i in 1:length(nonind))
                        if(sum(nonind==(nbr_test+1))==1){nbr_test=nbr_test+1}
                    }
                            
                    if(showresult) {cat("At least one alternative rejects H_k, the null hypothesis is rejected\n")}
                                
                }else{
                    T=0
                    if(showresult) {cat("The null hypothesis is accepted\n\n\n")}
                }
                                
            }# end while
            if(ktest==dim_X){k0=dim_X}else{k0=ktest}

            NBR=c(NBR,nbr_test)
            NBR_effect=c(NBR_effect,k0)

            relevant_variables=rbind(relevant_variables,c(ORDREBETA[1:nbr_test],rep(0,NBR[1]-nbr_test)))

    }#end alpha
        
    rownames(relevant_variables)=paste("alpha=",alpha)
    colnames(relevant_variables)=colnames(data)[relevant_variables[1,]]
    if(showresult){print("relevant variables:")
    print(relevant_variables)}
            
    quantile.all=array(0,c(length(alpha),maxqdep,dim_X)) #on y met tous les quantiles
    quantile.all[,,1:max(NBR_effect)]=quantile[,,1:max(NBR_effect)]
    rownames(quantile.all)=paste("alpha=",alpha)
    colnames(quantile.all)=paste("Hk,",0:(maxqdep-1))
    dimnames(quantile.all)[[3]]=paste("ktest=",(var_nonselect!=0):(dim_X-(var_nonselect==0)))

    quantile.all=quantile.all[,,1:max(NBR_effect),drop=FALSE]


    # fit a linear models on the relevant variables to get the coefficients
    Y.fitted=residuals=NULL
    coefficients=matrix(0,ncol=length(alpha),nrow=p)
    for(i in 1:length(alpha))
    {reg=lm(Y~data[,relevant_variables[i,]]-1)
        coefficients[relevant_variables[i,],i]=reg$coefficients
        reg$coefficients[-which(reg$coefficients!=0)]=0
        Y.fitted=cbind(Y.fitted,data[,relevant_variables[i,],drop=FALSE]%*%reg$coefficients)
        residuals=cbind(residuals,reg$residuals)
    }

    colnames(Y.fitted)=colnames(residuals)=paste("alpha=",alpha)
    rownames(coefficients)=colnames(data)
    colnames(coefficients)=alpha

    out=list(data=list(X=data,Y=Y,means.X=means.X,sigma.X=sigma.X,intercept=intercept),coefficients=coefficients,residuals=residuals,relevant_var=relevant_variables,fitted.values=Y.fitted,ordre=ORDRE,ordrebeta=ORDREBETA,kchap=NBR,quantile=quantile.all,call=match.call())

    out
    structure(out,class="mht")
}#fin procbol
