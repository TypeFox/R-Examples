#    <rank.btnr>
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

rank.btnr=function(data){

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
        gamma=numeric(c)
        change=1
        gm=gm=t(wm[,1:nmo])+wm[1:nmo,]
        wins=apply(wm[1:nmo,],1,sum)
        gind=(gm>0)
        squaregind=gind[,1:nmo]
        y=matrix(0,nmo,c)
        rus=y

        while(norm(as.matrix(change),"F")>1e-8){
              pi=exp(gamma)
              pius=matrix(rep(pi,c),ncol=c)
              piust=t(pius[,1:nmo])
              pius=pius[1:nmo,]
              rus[gind]=pius[gind]/(pius[gind]+piust[gind])
              y[gind]=gm[gind]*rus[gind]
              d2Q=-apply(y,1,sum)
              dL=wins+d2Q
              y[gind]=y[gind]*(1-rus[gind])
              d2L=-diag(apply(y,1,sum))
              temp=y[,-c]
              d2L[squaregind]=temp[squaregind]
            change=solve(d2L)%*%dL
            change=c(change,0)
            gamma=gamma-change
        }

        result=round(exp(gamma),4)
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


