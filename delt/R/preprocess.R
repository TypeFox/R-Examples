preprocess<-function(ssr,left,right,mean)
{
#ssr=excess mass

nodlkm<-length(ssr)
ssralip<-matrix(0,nodlkm,1)
#label<-matrix(1,nodlkm,1)

# Muodostetaan vektori, jossa kunkin solmun korkeus
kork<-matrix(0,nodlkm,1)
kork[1]<-1    #juuren korkeus on 1
i<-1
while (i<=nodlkm){
  if (right[i]>0){   #jos ei olla lehdessa
    kork[left[i]]<-kork[i]+1       #vasen lapsi yhta korkeammalla
    kork[right[i]]<-kork[i]+1      #oikea lapsi yhta korkeammalla
  }
  i<-i+1
}
tasolkm<-max(kork)     #tasojen lkm

k<-tasolkm             #aloitetaan korkeimmalta tasolta
while (k>0){
  i<-1
  while (i<=nodlkm){  #nodlkm
    if (kork[i]==k){
      if (right[i]==0){    #jos ollaan lehdessa
         ssralip[i]<-ssr[i] #alipuun ssr=solmun itsensa ssr
      }
      else if ((left[left[i]]==0) && (left[right[i]]==0)){
         minnu<-min(ssralip[left[i]],ssralip[right[i]])
         minnu<-min(minnu,ssr[i]) 
         if (ssralip[left[i]]==minnu){  #remove right
              ssralip[i]<-ssralip[left[i]]
              mean[right[i]]<-0
              left[right[i]]<-0
              right[right[i]]<-0
          }
          else 
             if (ssralip[right[i]]==minnu){  #remove left
                ssralip[i]<-ssralip[right[i]]
                mean[left[i]]<-0
                left[left[i]]<-0
                right[left[i]]<-0
             }
             else{
                left[i]<-0 
                right[i]<-0
                ssralip[i]<-ssr[i]
             }
      }  
      else{
         minnu<-min(ssralip[left[i]],ssralip[right[i]])
         minnu<-min(minnu,ssralip[left[i]]+ssralip[right[i]])
         minnu<-min(minnu,ssr[i])  
         if (minnu==ssralip[left[i]]+ssralip[right[i]]){
            ssralip[i]<-ssralip[left[i]]+ssralip[right[i]]
         }
         else
            if (ssralip[left[i]]==minnu){  #remove right
              ssralip[i]<-ssralip[left[i]]
              mean[right[i]]<-0
              left[right[i]]<-0
              right[right[i]]<-0
            }
            else 
              if (ssralip[right[i]]==minnu){  #remove left
                ssralip[i]<-ssralip[right[i]]
                mean[left[i]]<-0
                left[left[i]]<-0
                right[left[i]]<-0
              }
              else{
                left[i]<-0 
                right[i]<-0
                ssralip[i]<-ssr[i]
              }
      #
      }
    }
    i<-i+1
  }
  k<-k-1
}

return(list(right=right,left=left,S=ssralip,mean=mean))
}



