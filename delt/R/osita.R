osita<-function(n,wv,seed){
#
#if (wv>n) error
#
set.seed(seed)
#
koko<-floor(n/wv)
tulos<-matrix(0,koko+1,wv)
arpavec<-matrix(1,n,1)
i<-1
while (i<=n){
  arpavec[i]<-i
  i<-i+1
}
i<-1
while (i<=koko){
j<-1
while (j<=wv){
   uusidim<-n-((i-1)*wv+j)+1  
   arpa<-unidis(uusidim)
   tulos[i,j]<-arpavec[arpa]
   if (arpa==1){ 
       arpavec<-arpavec[2:uusidim]
   }
   else{
      if (arpa==uusidim){ 
        arpavec<-arpavec[1:(uusidim-1)]
      }
      else{ 
         arpavecnew<-matrix(0,uusidim-1,1)
         arpavecnew[1:(arpa-1)]<-arpavec[1:(arpa-1)]
         arpavecnew[arpa:(uusidim-1)]<-arpavec[(arpa+1):uusidim]
         arpavec<-arpavecnew
      }
   }
   j<-j+1
}
i<-i+1
} 
ylipitlkm<-n-wv*koko
j<-1
while (j<=ylipitlkm){
   uusidim<-n-(koko*wv+j)+1  
   arpa<-unidis(uusidim)
   tulos[koko+1,j]<-arpavec[arpa]
   if (arpa==1){ 
      arpavec=arpavec[2:uusidim]
   }
   else{
     if (arpa==uusidim){ 
        arpavec=arpavec[1:(uusidim-1)]
     }
     else{ 
         arpavecnew<-matrix(0,uusidim-1,1)
         arpavecnew[1:(arpa-1)]<-arpavec[1:(arpa-1)]
         arpavecnew[arpa:(uusidim-1)]<-arpavec[(arpa+1):uusidim]
         arpavec<-arpavecnew
     }
   } 
   j<-j+1
}
j<-ylipitlkm+1
while (j<=wv){
   tulos[koko+1,j]<-NA
   j<-j+1
}
#
return(tulos)
}





