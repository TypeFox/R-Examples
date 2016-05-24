# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction plot_pop_cluster_1d 
# fonction pour tracer les fréquences d'apparition 
# des mutations pour toutes les positions 
# des clusters positifs OU négatifs 
# pour toutes les séries d'évènement
#
# création : 09/04/14
# version du : 09/04/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.plot.cluster.all.1d<-function(x,events_series="all",...)
{ 
  cat("in progress ...\n")
  
  # -------------------------------------
  # vérifications et récupérer les données
  # ------------------------------------- 
  # longeur des séquences ou longeur des vecteurs
  # récupérer le facteur de normalisation (borne inf et sup de l'intervalle)
  # vérifier que toutes les bornes sont identiques
  if(is.null(x@param$n)) 
  {
    n<-x@param$normalisation_factor
    b_inf<-sapply(X=n,FUN=function(x) min(x))
    b_sup<-sapply(X=n,FUN=function(x) max(x))
    if(min(b_inf)!=max(b_inf) | min(b_sup)!=max(b_sup))
      stop("all events series must be on the same segment to plot with a vector of 'events_series'.",call.=F)
    n<-n[[1]]
  } else n<-c(0,x@param$n)
  
  # isoler les données à tracer
  d_total<-x@cluster$cluster_1d
  
  # éliminer séries d'évènement sans cluster détectés
  id<-!is.na(d_total[,1])
  d_total<-d_total[id,]
  
  
  # nombre de séries d'évènement
  events_series<-unique(d_total[,8])
  nb<-length(events_series)  
 

  # -------------------------------------
  # fréquence des évènements positions par
  # positions
  # -------------------------------------
  
  # déterminer si les positions possibles sont discrètes (sites) ou continues
  pos<-c(d_total[,1],d_total[,2])
  discrete<-all(((pos+1)/floor(pos+1))==1)
  
  # extraire les données en séparant cluster positifs et négatifs
  res_pos<-d_total[d_total[,5]==1,]
  res_neg<-d_total[d_total[,5]==0,]
  
  # vecteur des positions possibles de la séquence (nombre de sites)
  if(discrete)
  {
      position<-n[1]:n[2]    
  } else
  {
      position<-seq(from=n[1],to=n[2],length.out=2000)
  }
  nb_position<-length(position) 
  
  # fonction pour extraire en chaque position indices et nombres de clusters
  # x est une ligne du tableau res_pos ou res_neg
  extraire_position<-function(x)
  {
    if(discrete)
    {
	id<-x[1]:x[2] # ensemble des positions sur le vecteur concernées 
	              # les positions détectées sont comprises entre 0 et 1 donc il faut augmenter de +1
	indice[id+1]<<-indice[id+1]+x[3]    # incrémenter l'indice de cucala pour l'ensemble des positions
	frequence[id+1]<<-frequence[id+1]+1 # nombre d'indice pour pour l'ensemble des positions
    } else
    {
	id<-position>=x[1] & position<=x[2]
	id<-rank(position)[id]
	indice[id]<<-indice[id]+x[3]
	frequence[id]<<-frequence[id]+1
    }
    return(list(indice,frequence))
  }
    
  # appliquer la fonction 'extraire_position' sur 'res_pos' 
  indice<-rep(0,times=nb_position)    # vecteur indice de cucala en chaque positions
  frequence<-rep(0,times=nb_position) # vecteur fréquence des clusters en chaque positions  
  if(nrow(res_pos)>0) aa<-apply(res_pos,MARGIN=1,FUN=extraire_position)      
  denominateur<-frequence         # denominateur pour calculer la moyenne des indices en chaque position 
  denominateur[indice==0]<-1      # si aucun cluster en une position on divisera par 1 et non par zéro
  indice_pos<-indice/denominateur # indice_pos est une moyenne par position
  frequence_pos<-frequence
  position_pos<-position
  
  # appliquer la fonction 'extraire_position' sur 'res_neg' 
  indice<-rep(0,times=nb_position)    # vecteur indice de cucala en chaque positions
  frequence<-rep(0,times=nb_position) # vecteur fréquence des clusters en chaque positions  
  if(nrow(res_neg)>0)aa<-apply(res_neg,MARGIN=1,FUN=extraire_position)
  denominateur<-frequence         # denominateur pour calculer la moyenne des indices en chaque position 
  denominateur[indice==0]<-1      # si aucun cluster en une position on divisera par 1 et non par zéro
  indice_neg<-indice/denominateur # indice_neg est une moyenne par position
  frequence_neg<-frequence 
  position_neg<-position
  
  # éliminer les points chevauchants : plus de points en une position
  # que de séries d'évènements (nb)
  id<-frequence_pos<=nb
  position_pos<-position_pos[id]
  frequence_pos<-frequence_pos[id]
  indice_pos<-indice_pos[id]
  
  id<-frequence_neg<=nb
  position_neg<-position_neg[id]
  frequence_neg<-frequence_neg[id]
  indice_neg<-indice_neg[id] 
  
  # nombre final de positions
  nb_position_pos<-length(position_pos)
  nb_position_neg<-length(position_neg)
    
  # -------------------------------------
  # regroupement des fréquences
  # pour plus de 2000 fréquences
  # -------------------------------------
  if(nb_position_pos>2000)
  {
      id<-cut(position_pos,breaks=2000,include.lowest=T)
      frequence_pos<-tapply(frequence_pos,INDEX=id,FUN=max)
      position_pos<-tapply(position_pos,INDEX=id,FUN=min)
      indice_pos<-tapply(indice_pos,INDEX=id,FUN=mean)
      indice_pos[is.na(indice_pos)]<-0     
  }
  
  if(nb_position_neg>2000)
  {
      id<-cut(position_neg,breaks=2000,include.lowest=T)
      frequence_neg<-tapply(frequence_neg,INDEX=id,FUN=max)
      position_neg<-tapply(position_neg,INDEX=id,FUN=min)
      indice_neg<-tapply(indice_neg,INDEX=id,FUN=mean)
      indice_neg[is.na(indice_neg)]<-0     
  }  
  
  # nombre final de positions : remise à jour
  nb_position_pos<-length(position_pos)
  nb_position_neg<-length(position_neg)
  
  # -------------------------------------
  # éliminer NaN et infinis pour les indices
  # -------------------------------------
  indice_pos[is.na(indice_pos)]<-0
  indice_pos[is.infinite(indice_pos)]<-1.797693e+307 # 10 fois moins que le max 1.797693e+308
  indice_pos[indice_pos>=1.797693e+307]<-1.797693e+307
  
  indice_neg[is.na(indice_neg)]<-0
  indice_neg[is.infinite(indice_neg)]<-1.797693e+307   
  indice_neg[indice_neg>=1.797693e+307]<-1.797693e+307
  
  # -------------------------------------
  # couleurs des barres
  # -------------------------------------  
  palette<-colorRampPalette(colors=c("green","yellow","red"))(50)

  if(length(unique(indice_pos))>1)
  {
    # vérifier si le rapport entre min et max de l'indice est >10^2 (chercher une différence significative)  
    indice<-range(indice_pos[indice_pos>0])
    if(indice[2]/indice[1]<10^2) 
    {
	couleur_pos<-rep("green",times=nb_position_pos) # vecteur des couleurs
	log_pos<-NULL # ne pas tracer de légende
    }
    else
    {
	id<-indice_pos>0 # une couleur dans la palette que pour les positions avec cluster
	couleur_pos<-rep("white",times=nb_position_pos) # vecteur des couleurs
	# essayer avec une classification linéaire
	couleur<-as.character(cut(indice_pos[id],breaks=50,include.lowest=T,labels=palette))
	nb_classe<-length(table(couleur))
	log_pos<-0 # indicatrice pour échelle non log
	
	# si moins de 5 classes passer en classification log
	if(nb_classe<5)
	{
	  couleur<-as.character(cut(log(indice_pos[id]),breaks=50,include.lowest=T,labels=palette))
	  nb_classe<-length(table(couleur))
	  log_pos<-1 # indicatrice pour échelle log
	}
	couleur_pos[id]<-couleur
    }
  } 
  
  
  if(length(unique(indice_neg))>1)
  {
    # vérifier si le rapport entre min et max de l'indice est >10^2 (chercher une différence significative)  
    indice<-range(indice_neg[indice_neg>0])
    if(indice[2]/indice[1]<10^2) 
    {
	couleur_neg<-rep("green",times=nb_position_neg) # vecteur des couleurs
	exp_neg<-NULL # ne pas tracer de légende
    }
    else
    {
	id<-indice_neg>0 # une couleur dans la palette que pour les positions avec cluster
	couleur_neg<-rep("white",times=nb_position_neg) # vecteur des couleurs
	# essayer avec une classification linéaire
	couleur<-as.character(cut(indice_neg[id],breaks=50,include.lowest=T,labels=palette))
	nb_classe<-length(table(couleur))
	exp_neg<-0 # indicatrice pour échelle non exp
      
	# si moins de 5 classes passer en classification exp
	if(nb_classe<5)
	{
	  couleur<-as.character(cut(exp(indice_neg[id]),breaks=50,include.lowest=T,labels=palette))
	  nb_classe<-length(table(couleur))
	  exp_neg<-1 # indicatrice pour échelle exp
	}
	couleur_neg[id]<-couleur
    }
    
    
  } 
   
  # -------------------------------------
  # tracé des graphiques
  # -------------------------------------  
  # fonction pour tracer les barres
  tracer_barre<-function(x)
  {
    coordonnee<-matrix(c(x[1],0,x[1],x[2]),ncol=2,byrow=T)
    lines(coordonnee,col=x[3],lwd=2)
  }
  
  # --- tracé du graphique cluster positifs ---
  if((x@param$cluster_analysis=="both" | x@param$cluster_analysis=="positive") & length(unique(indice_pos))>1)
  {
      plot.new()
      plot.window(xlim=n,ylim=c(0,1.31))  
      apply(cbind(position_pos,frequence_pos/nb,couleur_pos),MARGIN=1,FUN=tracer_barre)
      
      # tracé des axes
      if(!is.null(x@param$n_sequences)) 
	  position<-round(quantile(n),0) else position<-quantile(n)
      axis(side=1,at=position)
      axis(side=2,at=seq(0,1,.25))

      # titres
      titre<-paste("events distribution, positives clusters",sep="")
      title(main=titre,xlab="positions",ylab="frequencies",font.main = 4)

      # ajout de la légende
      if(!is.null(log_pos))
      {
	    palette_legende<-colorRampPalette(colors=c("green","yellow","red"))(20)
	    if(log_pos==1) palette_legende<-colorRampPalette(colors=c("green","yellow","red"),bias=0.5)(20)
	    
	    xx1<-0.75*n[2]
	    xx2<-xx1+n[2]*0.08 
	    yy1<-1.05+(0:19)*0.015
	    yy2<-yy1+0.015
	    rect(xx1,yy1,xx2,yy2,col=palette_legende,border=NA)
	    rect(xx1,yy1[1],xx2,yy2[20],col=NA,border="black")

	    xx<-xx2+n[2]*0.0005
	    txt1<-format(min(indice_pos[indice_pos>0]),digits=3)
	    txt2<-format(max(indice_pos[indice_pos>0]),digits=3)
	    texte_legende<-c(txt1,txt2)
	    coord<-xy.coords(x=c(xx,xx),y=c(yy1[2],yy1[19])) # positions des labels
	    text(coord,labels=texte_legende,cex=1.1,pos=4)
      }
  }
   
  
  # --- tracé du graphique cluster négatifs ---
  if((x@param$cluster_analysis=="both" | x@param$cluster_analysis=="negative") &  length(unique(indice_neg))>1)
  {    
      dev.new()
      plot.new()
      plot.window(xlim=n,ylim=c(0,1.31))  
      apply(cbind(position_neg,frequence_neg/nb,couleur_neg),MARGIN=1,FUN=tracer_barre)
      
      # tracé des axes
      if(!is.null(x@param$n_sequences)) 
	  position<-round(quantile(n),0) else position<-quantile(n)
      axis(side=1,at=position)
      axis(side=2,at=seq(0,1,.25))

      # titres
      titre<-paste("events distribution, negatives clusters",sep="")
      title(main=titre,xlab="positions",ylab="frequencies",font.main = 4)

      # ajout de la légende
      if(!is.null(exp_neg))
      {
	  palette_legende<-colorRampPalette(colors=c("green","yellow","red"))(20)
	  if(exp_neg==1) palette_legende<-colorRampPalette(colors=c("green","yellow","red"),bias=2)(20)     
	  
	  xx1<-0.75*n[2]
	  xx2<-xx1+n[2]*0.08 
	  yy1<-1.05+(0:19)*0.015
	  yy2<-yy1+0.015
	  rect(xx1,yy1,xx2,yy2,col=palette_legende,border=NA)
	  rect(xx1,yy1[1],xx2,yy2[20],col=NA,border="black")

	  xx<-xx2+n[2]*0.0005
	  txt1<-format(-min(indice_neg[indice_neg>0]),digits=3)
	  txt2<-format(-max(indice_neg[indice_neg>0]),digits=3)
	  texte_legende<-c(txt1,txt2)
	  coord<-xy.coords(x=c(xx,xx),y=c(yy1[2],yy1[19])) # positions des labels
	  text(coord,labels=texte_legende,cex=1.1,pos=4)
      }
  }
  
}
