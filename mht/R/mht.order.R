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

#proc_ord=function(data,...){UseMethod("proc_ord")}
#proc_ord.default=function(data,Y,ordre,var_nonselect,alpha,IT,sigma,showresult,...)
mht.order=function(data,Y,ordre,var_nonselect,alpha,IT,sigma,showresult)
{
    #-----------------------------------
    #   data = Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.
    #   Y = Response variable of length n.
    #	ordre= variables rank. It can be a partial order. Defaut consider "data" has already order
    #	var_nonselect= Number of variables that are not submitted to feature selection; they are the first columns of the data matrix
    #	alpha= Type I error of the tests
    #	IT=Number of simulations to calculate the quantiles
    # 	sigma= positive number if the variance is known, otherwise 0. The statistics are different for known and unknown variance
    #	showresult=show the results of the multiple testing procedure
    #-----------------------------------------

    data = as.matrix(data)
    Y=as.numeric(Y)

    #-- validation des arguments --#
    if (length(dim(data)) != 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix.")

    # verifications de depart
    p=ncol(data)
    n=nrow(data)
    if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}
    if(missing(ordre)){ordre=1:p}
    if(missing(alpha)){alpha=c(0.1,0.05)}
    if(missing(showresult)){showresult=TRUE}
    if(missing(IT)){IT=10000}
    if(missing(sigma)){sigma=0}else{if (sigma<0) stop("sigma is a positive quantity")}

    if(p<3) stop("The function cannot run for less than 3 variables")

    ordreinit=ordre
    if(length(ordre)<p)
    {
        for(i in 1:p)
        {
            if(sum(i==ordre)==0){ordre=c(ordre,i)}
        } #on complete l'ordre par les variables restantes
    }
            
    #		-------------------------------------
    #			on scale la matrice de depart
    #		-------------------------------------

    #si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute
    temp=data.scale(data)
    data=temp$data
    intercept=temp$intercept
    means.X=temp$means.data
    sigma.X=temp$sigma.data
    p=ncol(data)
    
    if(intercept==FALSE)
    {
        ordre=c(1,ordre+1)#pour compenser l'ajout de l'intercept
    }
    
    dataa=data[,ordre]
    
    ordre=matrix(ordre,nrow=1)
    if(missing(var_nonselect))
    {
        var_nonselect=1
    }else{
        if(!intercept){var_nonselect=var_nonselect+1}
        if(var_nonselect<1){stop("var_nonselect has to be greater than 1")}
    }


    #		------------------------------------------
    #		  get an orthogonal family out of XI_ord
    #		------------------------------------------
    dec=decompbaseortho(dataa)# gives a orthogonal basis of  X(1),..X(p) where U[,i] is the decomposition of X(i) in the basis (Gram-Shmidt algorithm)
    nonind=dec$nonind #ex: nonind=4 means that the variable X(4) is a linear combination of X(1),X(2),X(3).
    U=dec$U     #U is the orthogonal basis of X(1),...,X(p).

    if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}
    dim_X=ncol(Uchap)			#nombre de variables utiles =dec$rank


    beta=lm(Y~Uchap-1)$coefficients # get the coefficients of the linear regression of Y onto the orthogonal basis U
    beta[-which(beta!=0)]=0


    #		----------------------------------------------------------------------------
    #			start of the sequential testing procedure
    #		----------------------------------------------------------------------------


    NBR=numeric(0) #record the number of hypotheses actually tested
    NBR_effect=numeric(0) #record the number of variables selected
    relevant_variables=numeric(0) #record the relevant variables

    quantile=matrix(0,length(alpha),min(n,dim_X)-1) #on y met tous les quantiles
    indice=0

    for(alph in 1:length(alpha)) #boucle sur alpha
    {

        ktest=var_nonselect-sum(nonind<=var_nonselect) # on commence par la premiere variable a selectionner
        
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
            if(ktest>indice)
            {
                quant=quantil_ord(n,dim_X,k=ktest,alpha,IT,sigma=sigma)#calcul du alpham
                quantile[,ktest]=quant$quantile
                indice=indice+1
            }
            
            alpha2=quantile[alph,ktest]

            bb=numeric(0)
            statistics.quantile=NULL
            maxq=log(min(n,dim_X)-ktest-1,2)
            for(m in 0:(maxq-1))
            {
                if(sigma==0)
                {
                    a=(n-ktest-2^m)*sum(beta[(ktest+1):(ktest+2^m)]^2)/(2^m*sum((Y-as.matrix(Uchap[,1:(ktest+2^m)])%*%beta[1:(ktest+2^m)])^2))
                    F=qf(1-alpha2,2^m,n-ktest-2^m)
                }else{
                    a=sum(beta[(ktest+1):(ktest+2^m)]^2)/sigma
                    F=qchisq(1-alpha2,2^m)
                }
                b=a>F #1 si on doit rejeter le test, 0 sinon
                bb=c(bb,b)
                statistics.quantile=cbind(statistics.quantile,c(a,F))
            }
            statistics.quantile=rbind(statistics.quantile,bb)
            rownames(statistics.quantile)=c("Statistics","Quantile","Rejected")
            colnames(statistics.quantile)=paste("S_(",ktest,",",0:(maxq-1),")",sep="")
            
            if(showresult){print(statistics.quantile)}

            if(sum(bb)>0)
            {
                T=1
                ktest=ktest+1
                
                nbr_test=nbr_test+1
                if(length(nonind)>0)
                {
                    for(i in 1:length(nonind))
                    if(sum(nonind==(nbr_test+1))==1){nbr_test=nbr_test+1}}
                
                if(showresult) {cat("At least one alternative rejects H_k, the null hypothesis is rejected\n")}
            }else{
                T=0
                if(showresult) {cat("The null hypothesis is accepted\n\n\n")}
            } #on rejete Hk s'il y a au moins un test qui rejete

        }
        if(ktest==dim_X){k0=dim_X}else{k0=ktest}


        NBR=c(NBR,nbr_test) #resultat contenant le nombre de variables selectionnees
        NBR_effect=c(NBR_effect,k0)

        relevant_variables=rbind(relevant_variables,c(ordre[1:nbr_test],rep(0,NBR[1]-nbr_test)))
    }#fin boucle sur alpha

    rownames(relevant_variables)=paste("alpha=",alpha)
    colnames(relevant_variables)=colnames(data)[relevant_variables[1,]]
    if(showresult){print("relevant variables:")
    print(relevant_variables)}

    if(length(alpha)==1){quantile.all=matrix(quantile[,1:max(NBR_effect)],nrow=1)}else{
    quantile.all=quantile[,1:max(NBR_effect),drop=FALSE]}

    rownames(quantile)=paste("alpha=",alpha)
    rownames(quantile.all)=paste("alpha=",alpha)
    colnames(quantile.all)=paste("ktest=",1:max(NBR_effect))	
    #fitted.values
    Y.fitted=residuals=NULL
    coefficients=matrix(0,ncol=length(alpha),nrow=p)
    for(i in 1:length(alpha))
    {
        reg=lm(Y~data[,relevant_variables[i,]]-1)
        coefficients[relevant_variables[i,],i]=reg$coefficients
        reg$coefficients[-which(reg$coefficients!=0)]=0
        Y.fitted=cbind(Y.fitted,data[,relevant_variables[i,],drop=FALSE]%*%reg$coefficients)
        residuals=cbind(residuals,reg$residuals)
    }

    colnames(Y.fitted)=colnames(residuals)=paste("alpha=",alpha)
    rownames(coefficients)=colnames(dataa)
    colnames(coefficients)=alpha
        
    out=list(data=list(X=data,Y=Y,means.X=means.X,sigma.X=sigma.X,intercept=intercept),coefficients=coefficients,residuals=residuals,relevant_var=relevant_variables,fitted.values=Y.fitted,ordre=ordreinit,ordrebeta=ordre,kchap=NBR,quantile=quantile.all,call=match.call())

    out
    structure(out,class="mht.order")
	
}