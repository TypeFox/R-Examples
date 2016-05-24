
# ' Zc colSums(Z)
missing_penalty<-function(nbclust_vect=nbclust_vect,Z=Z,M=M,n=n,p=p,Zc=Zc){
   penalty=0
      penalty_vect=nbclust_vect*3-1#3 parametres par classe mixmod mais somme des proportions a 1
   compl_vect=Zc+2#coef + constante+bruit
   used=rep(0,times=p)
   for (i in 1:n){
      used=0*used
      for(j in 1:p){
         if(M[i,j]==0){#variable observee
            if(Zc[j]==0){#variable a droite
               penalty=penalty+penalty_vect[j]
            }else{#variable a gauche
               penalty=penalty+compl_vect[j]#on compte la structure
            }
         }
      }
      quimankdroit=which(Z%*%(-(M[i,]-1))*M[i,]>0)
      penalty=penalty+penalty_vect[quimankdroit]
   }
   #variable a droite observee : mixmod compte a droite mais pas a gauche car sachant
   #variable a droite manquante : mixmod ne compte pas a droite car pas observee mais intervient dans la loi a gauche si la gauche est observee
   
   #vecteur binaire d'utilisation mis a jour a chaque ligne'
   penalty=log(n)*penalty/n#nombre moyen
   return(penalty)
}