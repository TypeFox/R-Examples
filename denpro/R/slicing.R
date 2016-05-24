slicing<-function(pcf,vecci,d1=1,d2=NULL)
# 1D slice: we calculate the slice parallel to the direction d1
# when d1=1, we calculate f(t,vecci) 
{
lenni<-length(pcf$value)
value<-matrix(0,lenni,1)
d<-length(pcf$N)
step<-matrix(0,d,1)
for (kk in 1:d) step[kk]<-(pcf$support[2*kk]-pcf$support[2*kk-1])/pcf$N[kk]

if  (is.null(d2)){  #1D slice

index<-matrix(0,lenni,1)

efek<-0
for (i in 1:lenni){

  currec<-matrix(0,2*d,1)
  for (kk in 1:d){
     currec[2*kk-1]<-pcf$support[2*kk-1]+pcf$down[i,kk]*step[kk]
     currec[2*kk]<-pcf$support[2*kk-1]+pcf$high[i,kk]*step[kk]
  }
  
  dimcal<-0
  onvalissa<-T
  j<-1
  while (j<=d){

     if (j!=d1){    
         ala<-currec[2*j-1]
         yla<-currec[2*j]
         if ((ala>vecci[j-dimcal]) || (yla<vecci[j-dimcal])) onvalissa<-F
     }
     else dimcal<-dimcal+1
     j<-j+1
  }
  if (onvalissa){
     efek<-efek+1
     value[efek]<-pcf$value[i]
     index[efek,1]<-pcf$high[i,d1]
  }
}

value<-value[1:efek]
index<-index[1:efek]
support<-pcf$support[(2*d1-1):(2*d1)]
N<-pcf$N[d1]
down<-matrix(0,efek,1)
high<-matrix(0,efek,1)
down[,1]<-index-1
high[,1]<-index
#down<-index-1
#high<-index

return(list(value=value,index=index,support=support,N=N,down=down,high=high))

}
else{ # 2D slice

if (is.null(d2)) d2<-2

down<-matrix(0,lenni,2)
high<-matrix(0,lenni,2)

efek<-0
for (i in 1:lenni){

  currec<-matrix(0,2*d,1)
  for (kk in 1:d){
     currec[2*kk-1]<-pcf$support[2*kk-1]+pcf$down[i,kk]*step[kk]
     currec[2*kk]<-pcf$support[2*kk-1]+pcf$high[i,kk]*step[kk]
  }
  
  dimcal<-0
  onvalissa<-T
  j<-1
  while (j<=d){

     if ((j!=d1) && (j!=d2)){    
         ala<-currec[2*j-1]
         yla<-currec[2*j]
         if ((ala>vecci[j-dimcal]) || (yla<vecci[j-dimcal])) onvalissa<-F
     }
     else dimcal<-dimcal+1
     j<-j+1
  }
  if (onvalissa){
     efek<-efek+1
     value[efek]<-pcf$value[i]
     down[efek,1]<-pcf$down[i,d1]
     down[efek,2]<-pcf$down[i,d2] 
     high[efek,1]<-pcf$high[i,d1]
     high[efek,2]<-pcf$high[i,d2] 
  }
}

value<-value[1:efek]
down<-down[1:efek,]
high<-high[1:efek,]
support<-c(pcf$support[(2*d1-1):(2*d1)],pcf$support[(2*d2-1):(2*d2)])

return(list(value=value,down=down,high=high,
support=support,N=c(pcf$N[d1],pcf$N[d2])))
} #else 2D slice

}




