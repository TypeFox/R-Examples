# ' MCMC
searchZmiss<-function(X=X,ini=NULL,maxit=10^5,M=NULL,plot=F,BIC_vrai=NULL,Z_vrai=NULL,mixmod=mixmod,nbclustmax=10,BIC_vide_vect=NULL, nett=T,mode=c("relax","rejet"),MH=T,p2max=Inf,rmax=5,nbit=1){
  #X est la matrice X sans la constante
  #Z est la matrice Z du formalisme priv?e de sa premi?re ligne (constante sur 1) et de sa premi?re colonne (identit?)
  #MH dit si on fait le saut vers le meilleur strict ou si on se contente de l'enregistrer
  mode=mode[1]
  X=1*as.matrix(X)
  p=ncol(X)
  n=nrow(X)
  if(is.null(M)){
    M=0*X
    M[is.na(M)]=1
    X[is.na(X)]=0
  }
  if(p2max<0){p2max=Inf}
  if(rmax<0){rmax=Inf}
  if(is.null(ini)){
    Z=matrix(0,ncol=p,nrow=p)
  }else{
    Z=matrix(ini)
  }
  Z_opt=Z
  if(is.null(mixmod)){
    mixmod=density_estimation(X=X,nbclustmax=nbclustmax,detailed=T)
    mixmod=mixmod_adapter(mixmod)
    BIC_vide_vect=mixmod$BIC_vect  
  }
  BIC_vide=sum(BIC_vide_vect)
  BIC_vect=BICZmiss(X=X,Z=Z_opt,Bic_vide_vect=BIC_vide_vect,intercept=T,mixmod=mixmod,nbit=nbit)
  BIC_opt=sum(BIC_vect)
  print(paste("BIC ini :",BIC_opt))
  
  #on initialise
  BIC=BIC_opt
  etape=0
  pas=0
  etape_max=round(maxit/(2*p-2))
  if(plot){
    compare=compare_struct(trueZ=Z_vrai,Zalgo=Z)
    courbe_BIC=rep(0,times=round(etape_max))
    courbe_BIC[1]=BIC
    courbe_BIC_opt=rep(0,times=round(etape_max))
    courbe_BIC_opt[1]=BIC_opt
    courbe_compl=rep(0,times=round(etape_max))
    courbe_compl[1]=sum(Z)
    courbe_nbbon1=rep(0,times=round(etape_max))
    courbe_nbbon1[1]=compare$nbbon1
    courbe_nbtrop=rep(0,times=round(etape_max))
    courbe_nbtrop[1]=compare$nbtrop
    courbe_nbmank=rep(0,times=round(etape_max))
    courbe_nbmank[1]=compare$nbmank
    courbe_bon_gauche=rep(0,times=round(etape_max))
    courbe_bon_gauche[1]=compare$bon_gauche
    courbe_faux_gauche=rep(0,times=round(etape_max))
    courbe_faux_gauche[1]=compare$faux_gauche
  }
  ptm<-proc.time()
  candidats=(1:p^2)[-seq(1,p^2,p+1)]
  
  nb_cand_loc=2*p-2
  BIC_loc_mat=matrix(0,ncol=3,nrow=(nb_cand_loc+1))
  meilleur_strict=F
  consecutif=T
  nbstat=0
  courbes=matrix()
  while(etape<maxit){
    ploc=sample(1:p,1)
    list_cand=rbind(cbind((1:p)[-ploc],ploc),cbind(ploc,(1:p)[-ploc])) 
    BIC_loc_mat[,c(2,3)]=0
    BIC_loc_mat[,1]=Inf
    BIC_loc_mat[1,]=c(BIC,0,0)
    if(mode=="rejet"){
      for (candidat in 1:nb_cand_loc){  
        Zloc=Z
        etape=etape+1
        i=list_cand[candidat,1]
        j=list_cand[candidat,2]
        Zloc[i,j]=1-Zloc[i,j] #on fait la modif
        if((Zloc[i,j] & sum(Zloc%*%Zloc)!=0)|(length(which(colSums(Zloc)!=0))>p2max)| (colSums(Zloc)[j]>rmax)){next}#croisement ou trop de ssreg ou ssreg trop compl, on passe                     
#         BIC_loc=sum(calcul_BIC2.2(Zloc=Zloc,X_appr=X,BIC_ini=BIC_vide_vect,BIC_Z=BIC_vect,Z=Z))
        BIC_loc=sum(BICZmiss(X=X,Z=Zloc,Bic_vide_vect=BIC_vide_vect,intercept=T,mixmod=mixmod,nbit=nbit,BicOld=BIC_vect,Zold=Z))
        BIC_loc_mat[candidat+1,]=c(BIC_loc,i,j)
#         if(is.na(BIC_loc)){print(paste("Sous-reg parfaite pour la variable",which(is.na(calcul_BIC2.2(Zloc=Zloc,X_appr=X,BIC_ini=BIC_vide_vect,BIC_Z=BIC_vect,Z=Z)))))}
      }
    }else if(mode=="relax"){
      for (candidat in 1:nb_cand_loc){  
        Zloc=Z
        etape=etape+1
        i=list_cand[candidat,1]
        j=list_cand[candidat,2]
        Zloc[i,j]=1-Zloc[i,j]
        if(as.numeric(Z[i,j])==1){#on supprime le croisement
          Zloc[,i]=0
          Zloc[j,]=0
        }  
        if((length(which(colSums(Zloc)!=0))>p2max)| (colSums(Zloc)[j]>rmax)){next}#si apr?s relax on est trop complexe, on passe
#         BIC_loc=sum(calcul_BIC2.2(Zloc=Zloc,X_appr=X,BIC_ini=BIC_vide_vect,BIC_Z=BIC_vect,Z=Z))
        BIC_loc=sum(BICZmiss(X=X,Z=Zloc,Bic_vide_vect=BIC_vide_vect,intercept=T,mixmod=mixmod,nbit=nbit,BicOld=BIC_vect,Zold=Z))
        BIC_loc_mat[candidat+1,]=c(BIC_loc,i,j)
        if(is.na(BIC_loc)){print(paste("Sous-reg parfaite pour la variable",which(is.na(BICZmiss(X=X,Z=Zloc,Bic_vide_vect=BIC_vide_vect,intercept=T,mixmod=mixmod,nbit=nbit,BicOld=BIC_vect,Zold=Z)))))}
      }
    }         
    #on termine l'?tape
    pas=pas+1      
    #on tire au sort entre stationnarit? et action
    loc=BIC_loc_mat[,1]
    if(min(loc[loc!=0])<BIC_opt){#on stocke, on pr?vient, et on y va
      #qui=which(BIC_loc_mat[1:length(loc),1]==min(BIC_loc_mat[1:length(loc),1]))[1]
      meilleur_strict=T
      Zloc=Z
      qui=which.min(loc)
      BIC_opt=BIC_loc_mat[qui,1]
      i=BIC_loc_mat[qui,2]
      j=BIC_loc_mat[qui,3]            
      if(mode=="rejet"){#donc on a pas eu de croisement
        Zloc[i,j]=1-Zloc[i,j] #on fait la modif      
      }else if(mode=="relax"){
        Zloc[i,j]=1-Zloc[i,j]
        if(as.numeric(Z[i,j])==1){#on supprime le croisement
          Zloc[,i]=0
          Zloc[j,]=0
        }                             
      }
      Z_opt=Zloc
      if(MH & meilleur_strict){#on saute
        meilleur_strict=F
        Z=Z_opt
        BIC=BIC_opt
        #print(paste("meilleur BIC",BIC))
      }
    }else{#si pas de meilleur, on tire au sort
      meilleur_strict=F#normalement inutile
      loc2=loc
      loc=loc[loc!=Inf]
      #on doit maintenant ?tablir la correspondance
      loc=loc-min(loc) #rend tout positif (deltabic) . Le min (ancien max) devient nul
      loc=exp(-0.5*loc)/sum(exp(-0.5*loc))
      loc=cumsum(loc)/sum(loc)#on cr?e les intervalles associ?s pour le tirage au sort
      station=which(loc>runif(1))[1]
      loc2[loc2!=Inf]=1:length(loc)
      station=which(loc2==station)#on va cherche l'indice correspondant      
      if(length(station)==0){
        print("aucun candidat r?alisable")
        station=1
      }
      if(station!=1){#alors pas de stationnarit?
        BIC=BIC_loc_mat[station,1]
        i=BIC_loc_mat[station,2]
        j=BIC_loc_mat[station,3]
        if(mode=="rejet"){#donc on a pas eu de croisement
          Z[i,j]=1-Z[i,j] #on fait la modif      
        }else if(mode=="relax"){
          Z[i,j]=1-Z[i,j]
          if(as.numeric(Z[i,j])==1){#on supprime le croisement               
            Z[,i]=0
            Z[j,]=0
          }                             
        }
      }else{nbstat=nbstat+1}  
    } 
    BIC_vect=BICZmiss(X=X,Z=Zloc,Bic_vide_vect=BIC_vide_vect,intercept=T,mixmod=mixmod,nbit=nbit)#completement sous-optimal

#     BIC_vect=calcul_BIC2.2(Zloc=Z,X_appr=X,BIC_ini=BIC_vide_vect)#voir si pas pr?f?rable de stocker l'etape precedente pour l'utiliser ici
    if(plot){
      compare=compare_struct(trueZ=Z_vrai,Zalgo=Z)
      courbe_BIC[pas]=BIC
      courbe_BIC_opt[pas]=BIC_opt
      courbe_compl[pas]=sum(Z)
      courbe_nbbon1[pas]=compare$nbbon1
      courbe_nbtrop[pas]=compare$nbtrop
      courbe_nbmank[pas]=compare$nbmank
      courbe_bon_gauche[pas]=compare$bon_gauche
      courbe_faux_gauche[pas]=compare$faux_gauche
    }
  } 
  print(proc.time()-ptm)
  if(nett){
       
  }
  if(plot){
    courbes=rbind(courbe_BIC,courbe_BIC_opt,courbe_compl,courbe_nbbon1,courbe_nbtrop,courbe_nbmank,courbe_bon_gauche,courbe_faux_gauche)
  }
  return(list(Z=Z_opt,BIC=BIC_opt,BIC_vide_vect=BIC_vide_vect,Zloc=Z,courbes=courbes,ptm=(proc.time()-ptm)[1],nbstat=nbstat))
}
