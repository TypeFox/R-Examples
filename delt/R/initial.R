initial<-function(ssr,left,right){
#Prunes a tree by removing useless splits.
#
#ssr, left, right are nodelkm-vectors
#tree is a list(ssr,left,right,...)
#
#Returnes left,right
#
#densplit saattaa muodostaa puun, jossa jokin alipuu
#on huonompi kuin ko. alipuun juuri. 
#Muokataan puuta siten, etta alipuut ovat vahintaan yhta hyvia.
#Algoritmi sama kuin Donoho 95, ts. 
#lahdetaan yhta tasoa lehtia ylempaa, edetaan taso kerrallan ylospain, 
#kustakin jaosta tarkistetaan, onko se hyodyllinen.
#Huom. oikea lapsi loydetaan hakemalla vasemman alipuun loppupiste
#ja siirtymalla yksi eteenpain: i:n oikea lapsi onn endpoint(tree,i+1)+1
#huom pyoristysvirhe vertailussa ??????????
#Kutsuu: poistamon, endpoint, sort.
#
ep<-0.0000001  #pyoristysvirheen huomioon ottaminen
#
nodlkm<-length(ssr)
#
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
# Muodostetaan vektori, jossa kustakin solmusta alkavan puun log-uskottavuus
#ssralip<-matrix(0,nodlkm,1)
#i<-nodlkm
#while (i>=1){
#  if (right[i]==0)    #jos ollaan lehdessa
#     ssralip[i]<-ssr[i] #alipuun ssr=solmun itsensa ssr
#  else ssralip[i]<-ssralip[left[i]]+ssralip[right[i]]
#            #muuten alipuun ssr=vasemman alipuun ssr + oik alip ssr
#  i<-i-1
#}
# Kaydaan puu lapi taso kerrallaan (miksi?) ja merkitaan muistiin poistettavat
#poistot<-matrix(0,nodlkm,1) #viite<-0   #poistojen maara
ssralip<-matrix(0,nodlkm,1)
tasolkm<-max(kork)     #tasojen lkm
k<-tasolkm             #aloitetaan korkeimmalta tasolta
while (k>0){
  i<-1
  while (i<=nodlkm){
    if (kork[i]==k){
      if (right[i]==0)    #jos ollaan lehdessa
         ssralip[i]<-ssr[i] #alipuun ssr=solmun itsensa ssr
      else{
          ssralip[i]<-ssralip[left[i]]+ssralip[right[i]]
            #muuten alipuun ssr=vasemman alipuun ssr + oik alip ssr
          if (ssralip[i]<=ssr[i]+ep){  
            #jos alipuu ei paranna, se poistetaan, huom pyoristysvirhe
             left[i]<-0 
             right[i]<-0
             #viite<-viite+1  #poistot[viite]<-i
        }
      }
    }
    i<-i+1
  }
  k<-k-1
}
#if (viite>=1){    #jos poistoja on tehty
#  poistot<-poistot[1:viite]
#}
#else{
#  poistot<-NULL
#}
return(list(right=right,left=left))
#return(list(dels=poistot,delnum=viite))
}












