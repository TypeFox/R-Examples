# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction plot_cluster_1d 
# fonction pour tracer les positions des clusters positifs et 
# négatifs pour une seule série d'évènement
#
# création : 17/12/13
# version du : 20/03/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.plot.cluster.one.1d<-function(x,events_series=1,...)
{
  # longeur de la séquence ou longeur du vecteur
  # récupérer le facteur de normalisation (borne inf et sup de l'intervalle)
  if(is.null(x@param$n)) 
    { n<-x@param$normalisation_factor[[events_series]] } else n<-c(0,x@param$n)
   
  # si une limite sur x est présente dans les paramètres suplémentaires
  argument<-list(...)
  if(length(argument)>0) if(match("xlim",names(argument))>0) n<-argument$xlim
  
  # isoler les données à tracer   
  d<-x@cluster$cluster_1d
  d<-d[d[,8]==events_series,]
  if(is.null(dim(d))) d<-matrix(d,nrow=1) # garder format matrix
   
   
  # nombre de clusters
  nb_cluster<-nrow(d)
  if(is.na(d[1,1])) stop(paste("no cluster plot for events serie ",
                              events_series,".",sep=""),call.=F)
  
  # -------------------------------------
  # tri sur les identifiants des clusters
  # -------------------------------------
  # traiter les clusters dans l'ordre de détection
  id<-sort(d[,6],index.return=T)
  d<-d[id$ix,]
  if(is.null(dim(d))) d<-matrix(d,nrow=1) # garder format matrix
   
  
  # -------------------------------------
  # mise à l'échelle des indices de Cucala
  # -------------------------------------
  indice_pos<-d[d[,5]==1,3]
  indice_neg<-d[d[,5]==0,3]
    
  if(length(indice_pos)==0) indice_pos<-1 else 
  {
      indice_pos<-log(indice_pos)    
  }
  if(length(indice_neg)==0) indice_neg<--1 else
  {
      indice_neg<--indice_neg/max(indice_neg)*mean(indice_pos)
  }

  # ajout d'une colonne à d pour y placer
  # les indices de Cucala mis à l'échelle
  d<-cbind(d,0)
  d[d[,5]==1,10]<-indice_pos
  d[d[,5]==0,10]<-indice_neg

  
  # -------------------------------------
  # calcul des proportions des longeurs
  # des clusters positifs et négatifs
  # -------------------------------------
  # facteur de normalisation
  normalisation<-n[2]-n[1] 
  pourcentage_longeur_pos<-sum(d[d[,5]==1,9])/normalisation*100
  pourcentage_longeur_neg<-sum(d[d[,5]==0,9])/normalisation*100
  pourcentage_longeur_pos<-round(pourcentage_longeur_pos,2)
  pourcentage_longeur_neg<-round(pourcentage_longeur_neg,2)
  
  # -------------------------------------
  # gestion des couleurs suivant 
  # les cluster positifs ou négatifs
  # -------------------------------------
  couleur<-rep("red",times=nrow(d))
  # cluster négatif
  couleur[d[,5]==0]<-"blue"           
  
  # -------------------------------------
  # tracé du graphique
  # -------------------------------------
  # ligne de base
  y0<-0
  
  plot.new()
  plot.window(xlim=n,ylim=range(d[,10]))
  rect(xleft=d[,1],xright=d[,2],ybottom=d[,10],
        ytop=y0,col=couleur)

  # tracé des axes
  if(!is.null(x@param$n_sequences)) 
      position<-round(quantile(n),0) else position<-quantile(n)
  axis(side=1,at=position)
  axis(side=2)
 
  # titres
  titre<-paste("cluster lengths: positives ",pourcentage_longeur_pos,"%, negatives ",pourcentage_longeur_neg,"%",sep="")
  
  if(x@param$cluster_analysis=="positive")
     titre<-paste("cluster lengths: positives ",pourcentage_longeur_pos,"%",sep="")
  if(x@param$cluster_analysis=="negative")
     titre<-paste("cluster lengths: negatives ",pourcentage_longeur_neg,"%",sep="") 
  
  sous_titre<-x@param$id[events_series]
  title(main=titre,sub=sous_titre,xlab="positions",ylab="Cucala index",font.main = 4)

}
