slicing.recs<-function(pa,vecci,d1=1,d2=2){

lenni<-length(pa$values)
d<-length(pa$recs[1,])/2

values<-matrix(0,lenni,1)
recs<-matrix(0,lenni,4)

efek<-0

for (i in 1:lenni){

  currec<-pa$recs[i,]
  
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
     values[efek]<-pa$values[i]
     recs[efek,1:2]<-pa$recs[i,(2*d1-1):(2*d1)]
     recs[efek,3:4]<-pa$recs[i,(2*d2-1):(2*d2)]

  }
}

values<-values[1:efek]
recs<-recs[1:efek,]

return(list(values=values,recs=recs))
}



