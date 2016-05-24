# ' initialization based on a wheight matrix (correlation or other)
# ' @param W Weight matrix with values between 0 and 1 ( e.g. : abs(cor(X)) )
# ' @param X the dataset (if BIC=TRUE)
# ' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
# '@param relax TRUE: constraint relaxation, FALSE: reject mode
# '@param p1max maximum complexity for a regression
# ' @param random if FALSE W is not used.
# ' @param mode candidates strategy
# '@param p2max maximum number of regressions 
# ' @param sorted boolean to sort the candidates based on W (best first)
# ' @param BIC boolean to use BIC for probabilities.
# ' @param nbclustmax parameter for calcul_BIC_mixmod if needed
# ' @export

Winitial<-function(W=W,X=NULL,p1max=NULL,BIC=F,Bic_null_vect=NULL,relax=T,random=T,nbclustmax=5,sorted=T,mode=c("all","sorted","multinom"),p2max=NULL){
  #W est la matrice des poids (entre 0 et 1)
  #X est la matrice des donnees (ne sert que si BIC=T)
  #Bic_null_vect est le vecteur bic vide (issu de mixmod par exemple)
  #relax indique si on fait la relaxation ou le rejet
  #p1max est le nombre max de 1 sur une colonne de Z
  #random indique si on fait le tirage au sort ou pas (si on utilise BIc on a moins besoin du tirage). sans tirage = sans tenir compte de W
  p=ncol(W)
  if(is.null(p2max)){p2max=p}
  diag(W)=0#on tue la diagonale au cas ou ce serait une matrice de corr?lation
  W=abs(W)
  mode=mode[1]
  Z=matrix(0,ncol=p,nrow=p)
  list_j=sample(p)#melange les entiers de 1 a p
  if(BIC==T &  !is.null(X)){
     if(is.null(Bic_null_vect)){
        Bic_null_vect=density_estimation(X=X)$BIC_vect 
     }
    BIC_opt_vect=Bic_null_vect
    BIC_opt=sum(Bic_null_vect)#initialisation de l'optimum
  }else{
    BIC=F#soit il l'etait deja soit il manque de quoi le calculer
  }
  lambda=0
  if(is.null(p1max)){p1max=p+1}#pour etre tranquille, revient a une borne infinie
  
  if(mode=="all"){#on teste tous les points
    for (j in 1:p){
      jloc=list_j[j]
      list_i=sample(p)#melange les entiers de 1 a p
      for (i in 1:p){
        if(list_i[i]!=jloc){#si on n'est pas sur la diagonale
          iloc=list_i[i]
          
          if(random){lambda=runif(1)}#on tire un nombre aleatoire entre 0 et 1
          if(!random | lambda<W[iloc,jloc] ){#si le poids est interessant ou bien pas d'alea
            if(sum(Z[,jloc])<p1max){#si on respecte la contrainte de complexite
              croisement=F
              if(sum(Z[jloc,])!=0 | sum(Z[,iloc])!=0){
                croisement=T
              }
              if(croisement==F | relax==T){#si pas croisement ou relax
                Zloc=Z
                Zloc[iloc,jloc]=1
                if(croisement==T){#donc croisement et relax
                  #faire la relax
                  Zloc[jloc,]=0
                  Zloc[,iloc]=0
                }
                if(length(which(colSums(Zloc)!=0))<=p2max){
                  if(BIC){#si on veut s'appuyer sur le bic
                    
                    BICloc_vect=BicZ(X=X,Z=Zloc,Bic_null_vect=Bic_null_vect,Bic_old=BIC_opt_vect,Zold=Z)
                    BICloc=sum(BICloc_vect)
                    if(BICloc<BIC_opt){#modification effective
                      Z=Zloc
                      BIC_opt_vect=BICloc_vect
                      BIC_opt=BICloc
                    }
                  }else{#on se fout du bic
                    Z=Zloc         
                  }
                }  
              }
            }
          }
        }
      }
    }
  }else if(mode=="sorted"){#on teste la colonne puis on y place un 1 selon les probas de la colonne (on priorise les meilleurs donc moins de bruit)
    for (j in 1:p){
      jloc=list_j[j]
      #on teste la colonne en fonction de sa plus grande proba
      if(random){
        go=rbinom(1,1,max(W[,jloc]))
      }else{
        go=1
      }
      if(go==1){#si la colonne est jugee interessante
        i=1
        list_i=1:p
        list_i=list_i[-jloc]#on enleve la diagonale
        list_i=list_i[order(W[-jloc,jloc],decreasing=T)]#on ordonne selon les poids
        while(go & length(list_i)>0){
          iloc=list_i[1]#on choisit le meilleur
          if(sum(Z[,jloc])<p1max){#si on respecte la contrainte de complexite
            croisement=F
            if(sum(Z[jloc,])!=0 | sum(Z[,iloc])!=0){
              croisement=T
            }
            if(croisement==F | relax==T){#si pas croisement ou relax
              Zloc=Z
              Zloc[iloc,jloc]=1
              if(croisement==T){#donc croisement et relax
                #faire la relax
                Zloc[jloc,]=0
                Zloc[,iloc]=0
              }
              if(length(which(colSums(Zloc)!=0))<=p2max){
                if(BIC){#si on veut s'appuyer sur le bic
                  
                  BICloc_vect=BicZ(X=X,Z=Zloc,Bic_null_vect=Bic_null_vect,Bic_old=BIC_opt_vect,Zold=Z)
                  BICloc=sum(BICloc_vect)
                  if(BICloc<BIC_opt){#modification effective
                    Z=Zloc
                    BIC_opt_vect=BICloc_vect
                    BIC_opt=BICloc
                  }
                }else{#on se fout du bic
                  Z=Zloc         
                }
              }
            }
            list_i=list_i[-1]
            if(length(list_i)>0){#s'il reste des candidats
              if(random){#si alea, on regarde si ca vaut la peine de continuer
                go=rbinom(1,1,max(W[list_i,jloc],0))#on regarde si on continue en fonction du meilleur restant     
              }else{#si pas d'alea on continue tant qu'il y a des candidats
                go=1
              }
            }else{
              go=0
            }
          }else{
            go=0
            list_i=c()#la contrainte de complexite est definitivement foutue donc on arrete
          }               
        } 
      }
    }  
  }else{#multinom
    for (j in 1:p){
      jloc=list_j[j]
      #on teste la colonne en fonction de sa plus grande proba
      if(random){
        go=rbinom(1,1,max(W[,jloc]))
      }else{
        go=1
      }
      if(go==1){#si la colonne est jugee interessante
        i=1
        list_i=1:p
        list_i=list_i[-jloc]#on enleve la diagonale
        while(go & length(list_i)>0){
          wloc=W[list_i,jloc]
          wloc=wloc/wloc
          qui=which(rmultinom(1,1,wloc)==1)
          iloc=list_i[qui]#on choisit selon une multinomiale ponderee par les poids normalises
          if(sum(Z[,jloc])<p1max){#si on respecte la contrainte de complexite
            croisement=F
            if(sum(Z[jloc,])!=0 | sum(Z[,iloc])!=0){
              croisement=T
            }
            if(croisement==F | relax==T){#si pas croisement ou relax
              Zloc=Z
              Zloc[iloc,jloc]=1
              if(croisement==T){#donc croisement et relax
                #faire la relax
                Zloc[jloc,]=0
                Zloc[,iloc]=0
              }
              if(length(which(colSums(Zloc)!=0))<=p2max){
                if(BIC){#si on veut s'appuyer sur le bic
                  
                  BICloc_vect=BicZ(X=X,Z=Zloc,Bic_null_vect=Bic_null_vect,Bic_old=BIC_opt_vect,Zold=Z)
                  BICloc=sum(BICloc_vect)
                  if(BICloc<BIC_opt){#modification effective
                    Z=Zloc
                    BIC_opt_vect=BICloc_vect
                    BIC_opt=BICloc
                  }
                }else{#on se fout du bic
                  Z=Zloc         
                }
              }
            }
            list_i=list_i[-qui]
            if(length(list_i)>0){#s'il reste des candidats
              if(random){#si al?a, on regarde si ca vaut la peine de continuer
                go=rbinom(1,1,max(W[list_i,jloc],0))#on regarde si on continue en fonction du meilleur restant     
              }else{#si pas d'alea on continue tant qu'il y a des candidats
                go=1
              }
            }else{
              go=0
            }
          }else{
            go=0
            list_i=c()#la contrainte de complexite est definitivement foutue donc on arrete
          }               
        }
        
      }
    }
    
  }
  return(Z)
}
