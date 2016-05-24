pruseqlev<-function(tree){

S<-tree$S
lu<-luo(tree)
p<-lu$p
G<-lu$G
g<-lu$g
N<-lu$N

left<-tree$left
right<-tree$right

alknodlkm<-length(tree$vec)  #solmujen lkm
alfalkm<-alknodlkm           #alfojen lkm <= lehtien lkm <= solmujen lkm

delnodebeg<-matrix(0,alfalkm,1)
delnodeend<-matrix(0,alfalkm,1)
delWbeg<-matrix(0,alfalkm,1)
delWend<-matrix(0,alfalkm,1)
leafs<-matrix(0,alfalkm,1)
alfa<-matrix(0,alfalkm,1)
loglik<-matrix(0,alfalkm,1)

delnodes<-matrix(0,alknodlkm,1)
delW<-matrix(0,alknodlkm,1)

alphamin<-0

if (N[1]==1){    #tree on triviaali 

  delnodes<-NULL
  delW<-NULL
  delnodebeg<-c(0)
  delnodeend<-c(0)
  delWbeg<-c(0)
  delWend<-c(0)
  leafs<-c(1)
  alfa<-c(0)
  loglik<-c(tree$ssr[1])
  tulos<-list(tree=tree,delnodes=delnodes,delW=delW,
              delnodebeg=delnodebeg,delnodeend=delnodeend,
              delWbeg=delWbeg,delWend=delWend,
              leafs=leafs,alfa=alfa,loglik=loglik)
}
else{ 
 k<-1           #tulee kertomaan alfojen lkm:n +1
 #remnodenum<-0  #poistettavien haarojen lkm
 #j<-0           #poistolaskuri

 delnodebeg[k]<-1
 delWbeg[k]<-1

 alpha<-alphamin

 while (N[1]>1){     #jos puu ei triviaali

   if (omaver(alpha,G[1])){    #(G[1]>alpha){
      leafs[k]<-N[1]
      alfa[k]<-alpha
      loglik[k]<-S[1] 
      alpha<-G[1]

      k<-k+1   #siirrytaan uuteen alfaan
      #j<-0     #poistolaskuri nollataan
 
      delnodebeg[k]<-delnodeend[k-1]+1 
      delnodeend[k]<-delnodeend[k-1]
      delWbeg[k]<-delWend[k-1]+1 
      delWend[k]<-delWend[k-1]


   }

   # Haetaan seuraavaksi solmu t, jonka alapuolelta voidaan katkaista.
   t<-1                        #aloitetaan juuresta
   while (omaver(G[t],g[t])){ #(G[t]<g[t]){  #jatk. kunnes g saav. miniminsa 
      if (omaver(G[left[t]],G[right[t]]))
                         #(omasam(G[t],G[left[t]])) 
           t<-left[t]                 #mennaan vasemmalle
      else t<-right[t]                #mennaan oikealle
   }

   delnodeend[k]<-delnodeend[k]+1
   #remnodenum<-remnodenum+1
   #j<-j+1
   delnodes[delnodeend[k]]<-t   

   # Tehdaan t:sta lehti
   N[t]<-1
   S[t]<-tree$ssr[t]
   G[t]<-NA                    #infty
   g[t]<-NA
   # Palataan juureen paivittaen N, S, g ja G.
   while (t>1){
      t<-p[t]
      ######################
      #S[t]<-S[left[t]]+S[right[t]]
         minnu<-min(S[left[t]],S[right[t]])
         minnu<-min(minnu,S[left[t]]+S[right[t]])
         minnu<-min(minnu,tree$ssr[t])  
         if (minnu==S[left[t]]+S[right[t]]){
            S[t]<-S[left[t]]+S[right[t]]
         }
         else
            if (S[left[t]]==minnu){  #remove right
              S[t]<-S[left[t]]
              delWend[k]<-delWend[k]+1
              delW[delWend[k]]<-right[t]
              #
              delnodeend[k]<-delnodeend[k]+1
              delnodes[delnodeend[k]]<-right[t] 
              #
              left[right[t]]<-0
              right[right[t]]<-0
              N[right[t]]<-1
              N[right[t]]<-1
              G[right[t]]<-NA
              G[right[t]]<-NA
              g[right[t]]<-NA
              g[right[t]]<-NA
            }
            else 
              if (S[right[t]]==minnu){  #remove left
                S[t]<-S[right[t]]
                delWend[k]<-delWend[k]+1
                delW[delWend[k]]<-left[t]
                #
                delnodeend[k]<-delnodeend[k]+1
                delnodes[delnodeend[k]]<-left[t] 
                #
                left[left[t]]<-0
                right[left[t]]<-0
                N[left[t]]<-1
                N[left[t]]<-1
                G[left[t]]<-NA
                G[left[t]]<-NA
                g[left[t]]<-NA
                g[left[t]]<-NA
              }
              else{
                S[t]<-tree$ssr[t]
                left[t]<-0 
                right[t]<-0
              }
      ###############################################
      if (left[t]==0){
                N[t]<-1
                G[t]<-NA
                g[t]<-NA
      }
      else{
         N[t]<-N[left[t]]+N[right[t]]
         g[t]<-(tree$ssr[t]-S[t])/(N[t]-1)
         G[t]<-omamindelt(g[t],omamindelt(G[left[t]],G[right[t]]))
      }

   }  #while t>1
 } #while N[t]>1

leafs[k]<-N[1]
alfa[k]<-alpha
loglik[k]<-S[1]  #palataan -log-likelista log-likeliin
alpha<-G[1]

leafs<-leafs[1:k]      #(k-1) kertoo alfojen maaran
alfa<-alfa[1:k]
loglik<-loglik[1:k]

delnodebeg<-delnodebeg[1:k]
delWbeg<-delWbeg[1:k]
delnodeend<-delnodeend[1:k]
delWend<-delWend[1:k]

delnodes<-delnodes[1:delnodeend[k]]
delW<-delW[1:delWend[k]]

tulos<-list(tree=tree,delnodes=delnodes,delW=delW,
            delnodebeg=delnodebeg,delnodeend=delnodeend,
            delWbeg=delWbeg,delWend=delWend,
            leafs=leafs,alfa=alfa,loglik=loglik)

}

return(tulos)
}







