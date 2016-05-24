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

dyadiqueordre=function(data,Y,m,maxordre,var_nonselect,showtest,showordre,random)
{
	
    #-----------------------------------	
    #   data = Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.
    #   Y = Response variable of length n.
    #	m= number of replications of the Lasso in the Bolasso technique. Only used when ordre is set to `bolasso'
    #	maxordre= Number of variables to be ordered in the first step of the procedure
    #	var_nonselect= Number of variables that are not submitted to feature selection; they are the first columns of the data matrix
    #	showordre=show the progress of the algorith
    #-----------------------------------------


    #		-------------------------------------
    #			checking entries
    #		-------------------------------------

    if(!is.matrix(data)) stop(" `data' has to be a matrix")
    if(is.matrix(Y)){
        if(ncol(Y)>1) stop("`Y' has to be a vector or a single column matrix")
    }
    
    Y=as.numeric(Y)
    n=nrow(data)
    p=ncol(data)
    if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}

    if(missing(m)){m=100}
    if(m<1) stop(" `m' has to be greater than 1")
    if(!missing(random))
    {
        if(!is.matrix(random)) stop("`random' has to be a matrix")
        if(ncol(random)!=m) stop("`random' has to have `m' columns")
    }

    if(missing(maxordre)){maxordre=min(n/2-1,p/2-1)}
    if(maxordre<0 | maxordre>min(n,p))
    {
        maxordre=min(n/2-1,p/2-1)
        message("maxordre should be in (0,p), it is set to ",maxordre)
    }

    if(missing(showordre)){showordre=TRUE}
    if(missing(showtest)){showtest=FALSE}

			
			#		-------------------------------------
			#			on scale la matrice de depart
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
            
    if(showordre){cat(colnames(data)[1]," ")}
    
    MU=numeric(0) #on va y mettre tous les mu que l'on teste
    COMPTEUR2=numeric(0) #on va y mettre tous les compteurs que l'on teste

    #on initialise le decoupage par le deuxieme mu d'un lasso classique (le premier etant tres grand)
    lasso1=glmnet(data,Y,alpha=1,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)))
        
    mu1=lasso1$lambda[2]
    muinf=0 #correspond au mu de gauche
    musup=mu1 #=mulasso #correspond au mu de droite


    if(var_nonselect>0){ordre=1:var_nonselect}else{ordre=numeric(0)}
    nonordonne=numeric(0)
    mu=musup
    mumil=musup #correspond au mu du milieu, celui qu'on teste pour savoir si on regarde a gauche, dans ce cas musup devient mumil, ou si on regarde a droite, dans ce cas muinf devient mumil.

    while(length(ordre)<maxordre)#on ordonne les variables jusqu'a en avoir maxordre
    {

    muinf=0
    musup=mu1
    mumil=musup
    longordo=0
    toto=0
    while(longordo!=1)
    {

        if(sum(MU==mu)==0)#si on a pas encore teste le mu
        {
            if(missing(random))
            {
                bol=bolasso(data,Y,mu,m,probaseuil=1,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)))
            }else{
                bol=bolasso(data,Y,mu,m,probaseuil=1,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)),random=random)
            }
        compteur2=bol$frequency

        MU=c(MU,mu)#on met tous les mu qu'on a deja teste
        COMPTEUR2=cbind(COMPTEUR2,compteur2)#on met les compteur2 de chaque mu teste
       
        if(showtest)
        {print(dim(COMPTEUR2))}
        
        }else{
            a=which(MU==mu) 			#si on a deja teste le mu
            compteur2=COMPTEUR2[,a]
        }	#on evite de refaire un calcul deja fait


        nonordonne=which(compteur2==1)#nous donne les indices des variables selectionnees

        #on regarde maintenant lesquelles etaient deja ordonnees
        ordonne=numeric(0)
        longnon=length(nonordonne)
        if(longnon>0)
        {
            for(i in 1:longnon)
            {
                if(sum(nonordonne[i]==ordre)==0)
                {
                    ordonne=c(ordonne,nonordonne[i])
                }#on recence les variables selectionnees qui n'ont pas encore ete ordonnees
            }
            longordo=length(ordonne) # nombre de difference entre ordre et nonordonne
        }else{
            longordo=0
        }
        

        #si il n'y a pas de nouvelles variables a ordonner, on diminue la penalite mu
        if(longordo==0)
        {
            musup=mumil
            mumil=(musup+muinf)/2
            mu=mumil
        }
            
        #si il y a plus qu'une unique variable a ordonner, on augmente la penalite mu pour n'en trouver qu'une seule
        if(longordo>1)
        {
            muinf=mumil
            mumil=(muinf+musup)/2
            mu=mumil
        }

        #on met une securite
        if((musup-muinf)<10e-10)
        {
            mu=muinf
            a=which(MU==mu) 			#si on a deja teste le mu
            compteur2=COMPTEUR2[,a] 	#on evite de refaire un calcul deja fait
            
            nonordonne=which(compteur2==1)#nous donne les indices des variables selectionnees
            
            #on regarde maintenant lesquelles etaient deja ordonnees
            ordonne=numeric(0)
            longnon=length(nonordonne)
            if(longnon>0)
            {
                for(i in 1:longnon)
                {
                    if(sum(nonordonne[i]==ordre)==0)
                    {
                        ordonne=c(ordonne,nonordonne[i])
                    }
                }
                longordo=length(ordonne) # nombre de difference entre ordre et nonordonne
            }else{
                longordo=0
            }
                
            break #nous sort du while s'il y a un probleme
        }

        if(showtest==TRUE) print(mu)

    }#fin while, on a donc une seule variable dans ordonne, c'est la variable suivante

    ordre=c(ordre,as.numeric(ordonne))

    if(showordre){cat(colnames(data)[as.numeric(ordonne)]," ")}

    }
    if(showordre) {cat("\n")}
    out=list(data=list(X=data,Y=Y,means.X=means.X,sigma.X=sigma.X),ordre=ordre,mu=MU,frequency=COMPTEUR2)
    out
    structure(out,class="bolasso")
}#fin function

