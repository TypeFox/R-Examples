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

quantilemht=function(data,k,alpha,IT,maxq,sigma)
{
    n=nrow(data)
    p=ncol(data)
        
    if(missing(alpha)){alpha=c(0.1,0.05)}
    if(missing(IT)){IT=1000}
    if(missing(maxq)){maxq=log(min(n,p)-1,2)}
    if(missing(sigma)){sigma=0}

    b=0
    nbrprob=0
    seuil=2
    while((sum(b==0)!=0)&(nbrprob<seuil)) # on recommence si a la fin on n'a pas assez de simulations pour obtenir une estimation du quantile
    {
        ijk=k
        if(ijk>0) #on doit construire une base de Vk
        {
            liste1=1:ijk
            dataa=as.matrix(data[,liste1]) 	#on doit construire une base ortho de dataa
            XI2dep=as.matrix(data[,-liste1])
                
            if(ijk>1)
            {
                dec=decompbaseortho(dataa)
                U=dec$U
            }else{U=svd(dataa)$u}
            U2dep=U
        }	

        sup=floor(maxq) #nombre max de tests que l'on fait
        #print(sup)
        F=matrix(0,IT*sup)
        if(k==0)
        {
            XI2dep=0
            U2dep=0
        }

        if(sigma>0)
        {
            epsilon=rnorm(n*IT,sd=sqrt(sigma))
        }else{
            epsilon=rnorm(n*IT)
        }

        a=quantiletest(k,data,XI2dep,U2dep,n,p,F,IT,maxq-1,sigma=sigma,e=epsilon)
        F=t(matrix(a$retour,sup,IT))
        #print(F[,sup])
        long=200
        lg_alpha=length(alpha)
        bV=array(0,c(long+1,sup,lg_alpha))
        aV=matrix(0,lg_alpha,sup)

        alphatest=matrix(0,long+1,lg_alpha)

        alphatest[1,]=alpha/log(n,2)
        for(j in 1:lg_alpha)
        {
            a=(alpha[j]-alpha[j]/(log(n,2)))/long
            for(i in 1:long)
            {
                alphatest[i+1,j]=alphatest[i,j]+a
            }
        }
        #alphatest

        FF=matrix(0,long+1,lg_alpha)
        for(m in 0:(maxq-1))
        {
            for(j in 1:lg_alpha)
            {
                for(i in 1:(long+1))
                {
                    FF[i,j]=quantile(F[,m+1],1-alphatest[i,j]) #quantile de U_{k,m}
                }
            bV[,m+1,j]=FF[,j] #chaque colonne=quantile, le numero de colonne est le m, la troisieme dim est le alpha
            }
        }

        FFF=matrix(0,long+1,lg_alpha)
        for(i in 1:IT)
        {
            FF=matrix(0,long+1,lg_alpha)
            for(m in 0:(maxq-1))
            {
                for(j in 1:lg_alpha)
                {
                    FF[,j]=FF[,j]+(F[i,m+1]>=bV[,m+1,j])	#=0 si la ieme iteration n'est pas rejete, >1 sinon
                }
            }

        #print(F)
        for(j in 1:lg_alpha)
        {
            FFF[,j]=FFF[,j]+(FF[,j]>0)}
        }

        b=numeric(0)
        for(j in 1:lg_alpha)
        {
            ind=which(FFF[,j]<=alpha[j]*IT)
            a=ind[length(ind)]
            #print(a)
            #b=c(b,a)
            if(length(a)==0){b=c(b,0)}else{b=c(b,a)}
        }
                    
        if(sum(b==0)!=0){nbrprob=nbrprob+1}#on incremente de 1 si on recommence

    }#fin while

    if(nbrprob==seuil){print(paste("number of simulations too low to estimate the ",alpha[which(b==0)],"-quantile, please increase IT"))}

    #break

    if(nbrprob<seuil)
    {
        for(j in 1:lg_alpha)
        {
            ind=which(FFF[,j]<=alpha[j]*IT)
            a=ind[length(ind)]
            aV[j,]=bV[a,,j]
        }
    }
    return(list(quantile=aV,nbrprob=nbrprob))
}#fin function
