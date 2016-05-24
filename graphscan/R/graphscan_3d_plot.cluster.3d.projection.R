# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction plot.cluster.3d 
# fonction pour tracer les clusters 3d avec une projection en 2d
# 
#
# création : 16/06/14
# version du : 16/06/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.plot.cluster.3d.projection<-function(x,indice="cucala",...)
{
    # choix de l'indice
    titre<-NULL
    if(indice=="cucala") titre<-"Cucala"
    if(indice=="kulldorff") titre<-"Kulldorff"
    if(is.null(titre))
      stop("argument 'indice' must be 'cucala' or 'kulldorff'.",call.=F)
    
    # position du SpatialPointsDataFrame dans la liste x@cluster en fonction de l'indice
    id<-grep(indice,names(x@cluster))
    if(length(id)==0)
      stop(paste("analysis with index '",indice,"' not present in graphscan object.",sep=""),call.=F)
      
    # récupérer les informations sur l'analysis et remise en forme
    info<-x@cluster$cluster_nd_description[id]    
    info<-gsub(pattern="- ",replacement="\n",x=info,fixed=T)
    
    # récupérer l'ensemble des cas 
    id_cas<-x@data$x$cases>0
    cas<-x@data$x[id_cas,]
    
    # rayon des cercles
    rayon<-x@cluster[[id]]$radius[1]
    
    # découper la fenêtre en 4 zones pour tracer
    par(mfrow=c(2,2))
    
    # --- graphe 1 y vs x ----
    # tracer tous les cas
    plot(cas@coords[,c(1,2)],pch=16,xlab="x",ylab="y")
        
    # tracer les cercles des clusters
    centre<-x@cluster[[id]]@coords[,c(1,2)]
    apply(X=centre,MARGIN=1,FUN=.cercle,rayon=rayon,couleur="lightgrey")
    
    # tracer les cas au centre des clusters
    points(x@cluster[[id]]@coords[,c(1,2)],pch=16,col="red")
      
    # affichage du titre
    title(main="y vs x")  
    
    # --- graphe 1 z vs x ----
    # tracer tous les cas
    plot(cas@coords[,c(1,3)],pch=16,xlab="x",ylab="z")
        
    # tracer les cercles des clusters
    centre<-x@cluster[[id]]@coords[,c(1,3)]
    apply(X=centre,MARGIN=1,FUN=.cercle,rayon=rayon,couleur="lightgrey")
    
    # tracer les cas au centre des clusters
    points(x@cluster[[id]]@coords[,c(1,3)],pch=16,col="red")
      
    # affichage du titre
    title(main="z vs x")  
    
    # --- graphe 1 z vs y ----
    # tracer tous les cas
    plot(cas@coords[,c(2,3)],pch=16,xlab="y",ylab="z")
        
    # tracer les cercles des clusters
    centre<-x@cluster[[id]]@coords[,c(2,3)]
    apply(X=centre,MARGIN=1,FUN=.cercle,rayon=rayon,couleur="lightgrey")
    
    # tracer les cas au centre des clusters
    points(x@cluster[[id]]@coords[,c(2,3)],pch=16,col="red")
      
    # affichage du titre
    title(main="z vs y")      
    
    # affichage des informations
    plot(0,axes=F,main="",xlab="",ylab="",col="white")
    text(x=0.57,y=0.5,pos=4,labels=info)
    
}