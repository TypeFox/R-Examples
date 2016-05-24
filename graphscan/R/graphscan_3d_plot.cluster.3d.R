# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction plot.cluster.3d 
# fonction pour tracer les clusters 3d 
# 
#
# création : 16/06/14
# version du : 16/06/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.plot.cluster.3d<-function(x,indice="cucala",sphere=TRUE,...)
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
    info<-gsub(pattern="\n",replacement="- ",x=info,fixed=T)
    
    # ouvrir le dispositif 3d
    open3d()
    bg3d(color="lightyellow")
    
    # tracer l'ensemble des cas 
    cas_id<-x@data$x@data$cases>0
    xyz<-x@data$x@coords[cas_id,]
    points3d(x=xyz[,1],y=xyz[,2],z=xyz[,3],col="black",size=5)

    # tracer centre des clusters
    xyz<-x@cluster[[id]]@coords
    rayon<-x@cluster[[id]]@data$radius[1]
    points3d(x=xyz[,1],y=xyz[,2],z=xyz[,3],col="red",size=10)
    # tracer sphères des clusters
    if(sphere==TRUE) spheres3d(x=xyz[,1],y=xyz[,2],z=xyz[,3],radius=rayon,col="grey",alpha=0.1,back="cull")
    
    # affichage des axes et des titres
    axes3d() 
    title3d(main=titre,sub=info,xlab="x",ylab="y",zlab="z")
}