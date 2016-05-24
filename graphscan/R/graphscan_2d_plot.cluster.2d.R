# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction plot.cluster.2d 
# fonction pour tracer les clusters 2d 
# 
#
# création : 24/03/14
# version du : 07/05/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.plot.cluster.2d<-function(x,map=NULL,indice="cucala",...)
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
      
    # récupérer les informations sur l'analysis
    info<-x@cluster$cluster_nd_description[id]
    
    # récupérer l'ensemble des cas 
    id_cas<-x@data$x$cases>0
    cas<-x@data$x[id_cas,]
    
    # tracer le contour map si nécessaire
    if(!is.null(map)) 
    { 
	if(class(map)!="SpatialPolygonsDataFrame" & class(map)!="SpatialPolygons")
	  stop("object 'map' must be of class 'SpatialPolygons' or 'SpatialPolygonsDataFrame'.",call.=F)
	plot(map) 
	# tracer tous les cas
	points(cas,pch=16)
    } else
    {   
	# tracer tous les cas
	plot(cas,pch=16)
    }
        
    # tracer les cercles des clusters
    rayon<-x@cluster[[id]]$radius[1]
    centre<-x@cluster[[id]]@coords
    apply(X=centre,MARGIN=1,FUN=.cercle,rayon=rayon,couleur="lightgrey")
    
    # tracer les cas au centre des clusters
    points(x@cluster[[id]],pch=16,col="red")
      
    # affichage titre
    title(main=titre,sub=info)  
}


# fonction pour tracer les clusters 2d
.cercle<-function(centre,rayon,couleur) 
{
  angle<-seq(0,2*pi,by=0.1)	
  xx<-rayon*cos(angle)+centre[1]
  yy<-rayon*sin(angle)+centre[2]   
  polygon(x=xx,y=yy,col=couleur,border=couleur)
}
