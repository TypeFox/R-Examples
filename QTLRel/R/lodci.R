
# search for LOD intervals #

.lodci<- function(llk,cv=0,lod=1.5,drop=3,index){
# llk: data.frame(dist,y,...)
# cv: a peak if llk$y>=cv
# drop: in case of multiple peaks, a peak if drop drop LOD on both sides
   nn<- nrow(llk)
   out<- NULL
   midx<- NULL
   nL<- nM<- nR<- 1
   xL<- xM<- xR<- llk$y[1]
   flag<- TRUE
   for(n in 1:nn){
      x0<- llk$y[n]
      if(x0>xM){
         nM<- n
         xM<- x0

         ii<- nL:n
         llkTmp<- llk[ii,]
         nL<- max(ii[llkTmp$y == min(llkTmp$y)])
         xL<- llk$y[nL]
         if(xM-xL >= drop && !flag){
            flag<- TRUE
         }
      }else{
         if(flag){
            if(xM-x0 >= drop && xM>= cv){
               ii<- nL:nM
               llkTmp<- llk[ii,]
               if(any(llkTmp$y <= xM-lod)){
                  nL<- max(ii[llkTmp$y <= xM-lod])
               }else nL<- max(ii[llkTmp$y == min(llkTmp$y)])
               xL<- llk$y[nL]

               ii<- nM:n
               llkTmp<- llk[ii,]
               if(any(llkTmp$y <= xM-lod)){
                  nR<- min(ii[llkTmp$y <= xM-lod])
               }else nR<- min(ii[llkTmp$y == min(llkTmp$y)])
               xR<- llk$y[nR]

               out<- rbind(out,c(llk$dist[nL],llk$dist[nR]))
               midx<- c(midx,index[nM])
               nL<- nM<- nR<- n
               xL<- xM<- xR<- x0
               flag<- FALSE
            }
         }else if(x0<=xL){
            nL<- n
            xL<- xM<- x0
         }
      }
   }
   if(flag && xM>=cv && xM-min(xL,llk$y[n]) >= drop){
      ii<- nL:nM
      llkTmp<- llk[ii,]
      if(any(llkTmp$y <= xM-lod)){
         nL<- max(ii[llkTmp$y <= xM-lod])
      }else nL<- max(ii[llkTmp$y == min(llkTmp$y)])
      xL<- llk$y[nL]

      ii<- nM:n
      llkTmp<- llk[ii,]
      if(any(llkTmp$y <= xM-lod)){
         nR<- min(ii[llkTmp$y <= xM-lod])
      }else nR<- min(ii[llkTmp$y == min(llkTmp$y)])
      xR<- llk$y[nR]

      out<- rbind(out,c(llk$dist[nL],llk$dist[nR]))
      midx<- c(midx,index[nM])
   }
   if(!is.null(out)) data.frame(lower=out[,1],upper=out[,2],index=midx) else out
}

lodci<- function(llk,cv=0,lod=1.5,drop=3){
   chrs<- unique(llk$ch)
   out<- data.frame(chr=NULL,lower=NULL,upper=NULL,index=NULL)
   for(n in 1:length(chrs)){
      idx<- c(1:nrow(llk))[llk$ch==chrs[n]]
      lk<- llk[idx,]
         lk<- data.frame(dist=lk$dist,y=lk$y)
      tt<- .lodci(lk,cv=cv,lod=lod,drop=drop,index=idx)
      if(!is.null(tt))
         out<- rbind(out,data.frame(chr=chrs[n],lower=tt$lower,upper=tt$upper,index=tt$index))
   }
   out
}


#############################
# old functions: not needed #
#############################

# local LOD support interval
flod<- function(llk,cv,lod=1.5,int=F){
   if(any(cv<lod)) stop("LOD too large.")
   if(max(llk$y)<cv) return(c(NA,NA))
   ii<- c(1:nrow(llk))[llk$y == max(llk$y)]
   ii<- ii[1]
   jj<- c(1:nrow(llk))[llk$y<llk$y[ii]-lod]
   cL<- cR<- NA
   if(any(jj<ii)){
      jL<- jj[jj<ii]
      jL<- max(jL)+1
      cL<- llk$dist[jL-1]
         if(int) cL<- cL + (llk$dist[jL]-llk$dist[jL-1])*(llk$y[ii]-lod-llk$y[jL-1])/(llk$y[jL]-llk$y[jL-1])
   }else cL<- llk$dist[1]
   if(any(jj>ii)){
      jR<- jj[jj>ii]
      jR<- min(jR)-1
      cR<- llk$dist[jR+1]
         if(int) cR<- cR + (llk$dist[jR]-llk$dist[jR+1])*(llk$y[ii]-lod-llk$y[jR+1])/(llk$y[jR]-llk$y[jR+1])
   }else cR<- llk$dist[nrow(llk)]

   c(cL,cR) #c(llk$dist[jL],llk$dist[jR])
}

fci<- function(llk,cv,lod){
   ci<- NULL
   dn<- c("lr","lrF2","lrF34")
   for(i in 1:length(cv)){
      llk$y<- llk[,dn[i]]
      for(j in 1:length(lod)){
         ci<- rbind(ci,flod(llk,cv[i],lod[j]))
      }
   }
   colnames(ci)<- c("lower","upper")
   dt<- rep(c("Integrated","F2","F34"),rep(length(lod),3))
   lod<- rep(lod,3)
   ci<- data.frame(data=dt,LOD=lod,ci)
   ci
}


