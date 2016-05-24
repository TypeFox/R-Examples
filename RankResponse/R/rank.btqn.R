#    <rank.btqn>
#    Copyright (C) <2014>  <Hsiuying Wang, Yu-Chun Lin>
#
#
#    This program is free software; you can redistribute it and/or modify
#
#    it under the terms of the GNU General Public License as published by
#
#    the Free Software Foundation; either version 2 of the License, or
#
#    (at your option) any later version.
#
# 
#
#    This program is distributed in the hope that it will be useful,
#
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#
#    GNU General Public License for more details.

rank.btqn=function(data){

        data=as.matrix(data)
        n=dim(data)[1]
        c=dim(data)[2]

        k=t(apply(data,1, function(x) x==1))
        B=matrix(0,c,c)
        wm=B;temp=B
        for(i in 1:n){
           B[k[i,],!k[i,]]=1
           wm=wm+B
           B=temp
        }

        nmo=c-1
        gamma=numeric(nmo)
        change=1
        gm=t(wm[,1:nmo])+wm[1:nmo,]
        wins=apply(wm[1:nmo,],1,sum)
        gind=(gm>0)
        y=matrix(0,nmo,c)
        z=y
        M=matrix(0,nmo,nmo)
        newchange=0
        i=0

        while(norm(as.matrix(change),"F")>1e-8){
              i=i+1
              pi=exp(gamma)
              pi[c]=1
                 pius=matrix(rep(pi,c),ncol=c)
                 piust=t(pius[,1:nmo])
                 pius=pius[1:nmo,]
                 y[gind]=gm[gind]/(pius[gind]+piust[gind])
                 z[gind]=pius[gind]*y[gind]
                 dL=wins-apply(z,1,sum)
                 newpi=wins/apply(y,1,sum)
              newgamma=log(newpi)
              oldchange=newchange
              newchange=newgamma-gamma
            if(i>1){
               g=olddL-dL
               s=-change+oldchange-newchange
               s_plus_Mg=s+M%*%g
               a=as.numeric(1/(t(g)%*%s_plus_Mg))
               M=M-a*(s_plus_Mg%*%t(s_plus_Mg))
            }
            olddL=dL
            change=newchange+M%*%dL
            gamma=gamma+change
        }

        result=round(c(exp(gamma),1),4)
        gamma=result/sum(result)
        names(gamma)=c(1:c)
        gamma_temp=as.numeric(names(sort(gamma,decreasing=T)))
        rank=numeric(c)
        
        for(i in 1:c)
        {
           rank[gamma_temp[i]]=i
        }
        result_rank=rbind(gamma,rank)
        return(result_rank)

}


