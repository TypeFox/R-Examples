pruseq<-function(tree){
#Forms a sequence of trees, which minimize likelihood-complexity
#criterion
#
#tree on list(val,vec,mean,nelem,ssr,left,right), ks densplit
#
#Result on list(tree,subtree,leafs,alfa,loglik)
#-tree on alkup. puu,
#-subtree on alfalkm*nodelkm-matriisi: osapuu annettu niitten solmujen 
#  indekseina, joista alkavat "tree":n alipuut eivat kuulu ko. osapuuhun.
#  nodelkm on yli alipuitten otettu maksimi poistettavien indeksien
#  maarasta
#-leafs,alfa,loglik are alfalkm-vectors: 
#    sis alipuuun lehtien lkm:n, alfan ja alipuun uskottavuuden

#alklehlkm<-leafnum(tree,1)     #lehtien lkm
alknodlkm<-length(tree$left)    #solmujen lkm
#nodelkm<-alknodlkm-alklehlkm
alfalkm<-alknodlkm             #alfojen lkm <= lehtien lkm <= solmujen lkm

delnodes<-matrix(0,alknodlkm,1)
delend<-matrix(0,alfalkm,1)      
leafs<-matrix(0,alfalkm,1)
alfa<-matrix(0,alfalkm,1)
loglik<-matrix(0,alfalkm,1)

# Initializing:

leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,alknodlkm,1)   #number of leaves in the tree whose root is i
R<--tree$ssr               #R(i) on noden i ssr eli -log likeli
p<-matrix(0,alknodlkm,1)   #parent
S<-matrix(0,alknodlkm,1)   #sen alipuun ssr -(log-likeli), jonka juuri i
g<-matrix(0,alknodlkm,1)   #(R(i)-S(i))/(N(i)-1), R(i) on noden i ssr
G<-matrix(0,alknodlkm,1)   #min{g(t),G(l(t)),G(r(t))}
t<-alknodlkm
while (t>=1){
  if (!is.na(leafloc[t]) && (leafloc[t]==1)){  #l(t)=0 eli ollaan lehdessa 
   N[t]<-1
   S[t]<--tree$ssr[t]                        
   G[t]<-NA                 #\infty
  }
  else  if (!is.na(leafloc[t])){
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     S[t]<-S[tree$left[t]]+S[tree$right[t]]  #S[t+1]+S[endpoint(tree,t+1)+1]    
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
     g[t]<-(R[t]-S[t])/(N[t]-1)
     G[t]<-omamindelt(g[t],omamindelt(G[tree$left[t]],G[tree$right[t]]))
  }
  t<-t-1
}

alphamin<-0

if (N[1]==1){    #tree on triviaali 

  delnodes<-NULL
  delend<-c(0)
  leafs<-c(1)
  alfa<-c(0)
  tulos<-list(tree=tree,delnodes=delnodes,delend=delend,leafs=leafs,alfa=alfa,
              loglik=loglik)
}
else{ 
 k<-1  #tulee kertomaan alfojen lkm:n +1
 delend[k]<-0
 j<-0  #poistolaskuri
 remnodenum<-0    #poistettavien haarojen lkm
 alpha<-alphamin
 while (N[1]>1){     #jos puu ei triviaali
   if (omaver(alpha,G[1])){    #(G[1]>alpha){
      leafs[k]<-N[1]
      alfa[k]<-alpha
      loglik[k]<--S[1]  #palataan -log-likelista log-likeliin
      alpha<-G[1]
      if (k>=2) delend[k]<-delend[k-1]+j
      k<-k+1   #siirrytaan uuteen alfaan
      j<-0     #poistolaskuri nollataan
   }
   #Haetaan seuraavaksi solmu t, jonka alapuolelta voidaan katkaista.
   t<-1                        #aloitetaan juuresta
   while (omaver(G[t],g[t])){ #(G[t]<g[t]){  #jatk. kunnes g saav. miniminsa 
      if (omaver(G[tree$left[t]],G[tree$right[t]]))
                         #(omasam(G[t],G[tree$left[t]])) 
         t<-tree$left[t]                 #mennaan vasemmalle
      else t<-tree$right[t]              #mennaan oikealle
   }
   remnodenum<-remnodenum+1
   delnodes[remnodenum]<-t   
   j<-j+1
   #Tehdaan t:sta lehti
   N[t]<-1
   S[t]<-R[t]
   G[t]<-NA                    #infty
   #Palataan juureen paivittaen N, S, g ja G.
   while (t>1){
      t<-p[t]
      S[t]<-S[tree$left[t]]+S[tree$right[t]]
      N[t]<-N[tree$left[t]]+N[tree$right[t]]
      g[t]<-(R[t]-S[t])/(N[t]-1)
      G[t]<-omamindelt(g[t],omamindelt(G[tree$left[t]],G[tree$right[t]]))
   }
 }

leafs[k]<-N[1]
alfa[k]<-alpha
loglik[k]<--S[1]  #palataan -log-likelista log-likeliin
alpha<-G[1]
if (k>=2) delend[k]<-delend[k-1]+j
k<-k+1

leafs<-leafs[1:(k-1)]      #(k-1) kertoo alfojen maaran
alfa<-alfa[1:(k-1)]
loglik<-loglik[1:(k-1)]
delnodes<-delnodes[1:remnodenum]
delend<-delend[1:(k-1)]
tulos<-list(tree=tree,delnodes=delnodes,delend=delend,leafs=leafs,alfa=alfa,
loglik=loglik)
}

return(tulos)
}







