# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction plot_pop_cluster_1d 
# fonction pour tracer les positions des clusters positifs et 
# négatifs pour un ensemble de séries d'évènement
#
# création : 17/12/13
# version du : 09/04/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.plot.cluster.serie.1d<-function(x,events_series=1:6,...)
{ 
  
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
 
  # nombre de séries d'évènement
  nb<-length(events_series)    
  
  # données totales extraites
  d_total<-x@cluster$cluster_1d
  
  # isoler les données à tracer
  # filtrer sur events_series
  d<-NULL
  for(i in events_series) d<-rbind(d,d_total[d_total[,8]==i,])

  # -------------------------------------
  # calcul position en y des différentes séries
  # ------------------------------------- 
  k<-rep((1:nb),time=table(d[,8]))
  d<-cbind(d,k+(d[,5]==1)*0.05)
  
  # -------------------------------------
  # création des identifants des différentes séries
  # -------------------------------------   
  id<-x@param$id[events_series]
  id<-rep(id,time=table(d[,8]))
  id[is.na(d[,1])]<-paste(id[is.na(d[,1])]," (no cluster)",sep="")
  id<-unique(id)
  
  # -------------------------------------
  # fonction pour tracer les clusters
  # -------------------------------------
  # palette des couleurs suivant les cluster positifs ou négatifs
  palette_cluster<-c("blue","red")
  
  tracer_cluster<-function(x)
  {
    coordonnee<-matrix(x[c(1,10,2,10)],ncol=2,byrow=T)
    lines(coordonnee,col=palette_cluster[x[5]+1],lwd=2)
  }
  
  # -------------------------------------
  # tracé du graphique
  # -------------------------------------
  plot.new()
  plot.window(xlim=n,ylim=c(0,nb+1))
  apply(d[!is.na(d[,1]),],MARGIN=1,FUN=tracer_cluster)

  # tracé des axes
  if(!is.null(x@param$n_sequences)) 
      position<-round(quantile(n),0) else position<-quantile(n)
  axis(side=1,at=position)
  
  # labels des identifant
  identifiant_x<-min(n)
  identifiant_y<-(1:nb)+0.3
  text(x=identifiant_x,y=identifiant_y,labels=id,pos=4)

  # titres
  titre<-paste("cluster localisations",sep="")
  title(main=titre,xlab="positions",ylab="",font.main = 4)

}
