# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction barplot_cluster_length_1d 
# fonction pour tracer la distribution des longeurs des clusters
# création : 17/12/13
# version du : 05/09/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

setMethod(f="barplot",signature="graphscan",definition=function(height,...)
{
    if(!is.null(height@cluster$cluster_1d_description))
    {
      # extraire les données de l'objet graphscan
      d<-height@cluster$cluster_1d_description[,c(3,4)]
      
      # couleur des barres
      couleur_pos<-colorRampPalette(colors=c("red1","red4"))(12)
      couleur_neg<-colorRampPalette(colors=c("blue1","blue4"))(12)
      par(mfrow=c(1,2))
      
      # cluster positifs
      a<-hist(d[,1],freq=T,main="",xlab="positive",ylab="number of events series",breaks=10,col=couleur_pos,axes=F)
      # positions et étiquettes de l'axe horizontal 
      position<-seq(from=min(a$breaks),to=max(a$breaks),length.out=5)
      etiquette<-paste(round(position,2),"%",sep="")
      axis(side=1,at=position,labels=etiquette)
      axis(side=2)

      # affichage du titre
      title(main="Cluster lengths",font.main=4,cex.main=0.9)
      
      # cluster négatifs
      a<-hist(d[,2],freq=T,main="",xlab="negative",ylab="number of events series",breaks=10,col=couleur_neg,axes=F)
      position<-seq(from=min(a$breaks),to=max(a$breaks),length.out=5)
      etiquette<-paste(round(position,2),"%",sep="")
      axis(side=1,at=position,labels=etiquette)
      axis(side=2)
      
      # affichage du titre
      titre<-paste(height@param$n_events_series," events series.",sep="")
      title(main=titre,font.main=4,cex.main=0.9)
      
    } else stop("graphscan object no containt cluster.",call.=F)
})



