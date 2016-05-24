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

refit.mht.order=function(object,Ynew,ordrenew,IT,var_nonselect,sigma,showresult,...)
{


	if(missing(Ynew)){stop('Ynew is missing')}
    if(missing(showresult)){showresult=TRUE}
	if(missing(sigma)){sigma=0}
	if(missing(IT)){IT=10000}
    
    
    #		-----------------------------------------
    #			recover some parameters from object
    #		-----------------------------------------
    
	alpha=as.numeric(colnames(object$coefficients))
	data=object$data$X
	Y=object$data$Y	
    means.X=object$data$means.X
    sigma.X=object$data$sigma.X
    intercept=object$data$intercept
    
    n=nrow(object$data$X)
    p=ncol(object$data$X)
    
	if(length(Ynew)!=n){stop(" 'data' and 'Ynew' must have the same length ")}



    if(missing(ordrenew))
    {
        ordrenew=object$ordre
    }
    if(intercept==FALSE)
    {
        ordrenew=c(1,ordrenew+1)#pour compenser l'ajout de l'intercept
    }
    

    if(length(ordrenew)<p)
    {
        for(i in 1:p)
        {
            if(sum(i==ordrenew)==0){ordrenew=c(ordrenew,i)} #on complete l'ordre par les variables restantes
        }
    }
    dataa=data[,ordrenew]

        
    ORDREBETA=ordrenew
    
    Y=Ynew
    ORDREBETA2=matrix(object$ordrebeta,ncol=p)
    
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


    quantile.all=array(0,c(length(alpha),dim_X-1,nrow(object$ordrebeta)+1))
    quantile.all[,1:dim(object$quantile)[2],1:nrow(object$ordrebeta)]=object$quantile
    #print(I)
    INDICE=matrix(0,nrow=2,ncol=nrow(object$ordrebeta))
    #pour chaque ordre i on a INDICE[1,i] qui nous donne var_nonselect de l'ordre i et INDICE[2,i] qui nous donne le ktest max
    for(i in 1:nrow(object$ordrebeta))
    {
        j=1
        while(sum(quantile.all[,j,i]^2)==0) {j=j+1}
        INDICE[1,i]=j-1
        while(sum(quantile.all[,j,i]^2)!=0) {j=j+1
        if(j>ncol(object$quantile)){break}}
        INDICE[2,i]=j-1			
    }
                              
    if(missing(var_nonselect)) 
    {
        var_nonselect=INDICE[1,i]+1
    }else{	
        if(var_nonselect<1){stop("var_nonselect has to be greater than 1")}
    }


    #		----------------------------------------------------------------------------
    #			on test pour separer les variables en deux: les bonnes et les autres
    #		----------------------------------------------------------------------------
    
    NBR=numeric(0) #record the number of hypotheses actually tested
    NBR_effect=numeric(0) #record the number of variables selected
    relevant_variables=numeric(0) #record the relevant variables

    quantile=matrix(0,length(alpha),min(n,dim_X)-1) #on y met tous les quantiles
    indice=0
    calcul=numeric(0)
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
                #on a indice de longueur compteur-1, quantile_compteur avec les quantiles, var_select de dim compteur-1, et ORDREBETA2
            if(ktest>indice)
            {
                abc=nrow(ORDREBETA2)#dim de indice/ORDREBETA2/var_select
                i=0
                I=numeric(0)
                TT=0
                while((TT==0)&&(i<abc))#dim(ORDREBETA2)=abc*maxordre
                {
                    i=i+1
                    a=numeric(0)
                    K=numeric(0)
                    for(j in 1:nbr_test)
                    {
                        a=c(a,sum(ORDREBETA[j]==ORDREBETA2[i,1:nbr_test]))
                    }
                    if((sum(a)==nbr_test)&&(INDICE[2,i]>=ktest)&&(INDICE[1,i])<ktest)
                    {
                        TT=1
                        #K=c(K,ktest)
                        I=c(I,i)
                    }
                }

                if(length(I)!=0)
                {
                    quantile[,ktest]=quantile.all[,ktest,I[length(I)]]#object$quantile[,,ktest]#get(paste("quantile",I,sep="_"))[,,ktest]
                    calcul=c(calcul,0)}else{
                    quant=quantil_ord(n,dim_X,k=ktest,alpha,IT,sigma=sigma)#calcul du alpham
                    quantile[,ktest]=quant$quantile
                    calcul=c(calcul,1)
                }
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
                if(showresult) {cat("The null hypothesis is accepted\n\n\n") }
            } #on rejete Hk s'il y a au moins un test qui rejete
                                
        }
        if(ktest==dim_X){k0=dim_X}else{k0=ktest}
                                    
        NBR=c(NBR,nbr_test) #resultat contenant le nombre de variables selectionnees
        NBR_effect=c(NBR_effect,k0)		
        relevant_variables=rbind(relevant_variables,c(ORDREBETA[1:nbr_test],rep(0,NBR[1]-nbr_test)))

    }#fin boucle sur alpha
        
    rownames(relevant_variables)=paste("alpha=",alpha)
    colnames(relevant_variables)=colnames(data)[relevant_variables[1,]]
    if(showresult){print("relevant variables:")
        print(relevant_variables)}


    if(sum(calcul)!=0)
    {
        quantile.all[,,nrow(object$ordrebeta)+1]=quantile

        if(length(alpha)==1)
        {
            quantile.all=array(quantile.all[,1:max(NBR_effect,dim(object$quantile)[2]),],c(1,max(NBR_effect,dim(object$quantile)[2]),nrow(object$ordrebeta)+1))
        }else{
            quantile.all=quantile.all[,1:max(NBR_effect,dim(object$quantile)[2]),]
        }

                    
        ORDREBETA=rbind(ORDREBETA2,ORDREBETA)

        rownames(quantile.all)=paste("alpha=",alpha)
        colnames(quantile.all)=paste("ktest=",1:max(NBR_effect,dim(object$quantile)[2]))	
        dimnames(quantile.all)[[3]]=paste("ordrebeta=",1:dim(quantile.all)[3])
    }else{
        quantile.all=object$quantile
        ORDREBETA=object$ordrebeta
    }

    #fitted.values


    Y.fitted=residuals=NULL
    coefficients=matrix(0,nrow=p,ncol=length(alpha))
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

    out=list(data=list(X=data,Y=Y,means.X=means.X,sigma.X=sigma.X,intercept=intercept),coefficients=coefficients,residuals=residuals,
    relevant_var=relevant_variables,fitted.values=Y.fitted,ordre=object$ordre,
    ordrebeta=ORDREBETA,kchap=NBR,quantile=quantile.all,call=match.call(),call.old=object$call)

    out
    structure(out,class="mht.order")

		
}#fin refit proc_ord
