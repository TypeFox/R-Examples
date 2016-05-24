# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction graphscan_plot 
# fonction générale pour tracer clusters 1d et nd
#
# création : 18/03/14
# version du : 07/07/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

# --------------------------------
# fonction .plot.graphscan
# --------------------------------
.plot.graphscan<-function(x,events_series=1,map=NULL,indice="cucala",sphere=TRUE,projection=FALSE,...)
{
  # --------------------------------------------
  # vérifications générales sur x (objet graphscan)
  # --------------------------------------------
  # 1) vérifier que la recherche de clusters a bien été effectuée
  if(x@param$dimension=="1d")
  {
    if(is.null(x@cluster$cluster_1d)) 
      stop("cluster analysis must be done.",call.=F)
  }
    
  if(x@param$dimension=="2d" | x@param$dimension=="3d")
  {
     if(length(x@cluster)==0)
       stop("cluster analysis must be done.",call.=F)
  } 
    
  # 2) vérifier qu'au moins un cluster a bien été détecté
  if(x@param$dimension=="1d")
  {
    if(class(x@cluster$cluster_1d_description)!="matrix")
      stop("no cluster has been detected in the data.",call.=F)
  }
   
   if(x@param$dimension=="2d" | x@param$dimension=="3d")
  {
     if(is.null(x@cluster$cluster_nd_cucala) & is.null(x@cluster$cluster_nd_kulldorff))
       stop("no cluster has been detected in the data.",call.=F)
  }
  
  # 3) vérification sur le paramètre events_series
  # rendre unique events_series
  events_series<-unique(events_series)
  
  # vérifier l'absence de NA dans events_series
  if(any(is.na(events_series))) 
     stop("'NA' present in 'events_series'.",call.=F)
  
  # --------------------------------------------
  # graphe 1d
  # --------------------------------------------  
  if(x@param$dimension=="1d")
  {    
    # lancer fonction graphe 1d
    if(length(events_series)==1 & min(events_series)!="all") # une seule série d'évènement
    {
      # vérifications sur events_series
      if(is.numeric(events_series) & events_series>x@param$n_events_series)
      stop(paste("argument 'events_series' must be less than or equal to ",x@param$n_events_series,sep=""),call.=FALSE)

      if(is.numeric(events_series) & events_series<1)
      stop("argument 'events_series' must be greater than 0.",call.=F) 
      
      if(is.character(events_series) & is.na(match(events_series,x@param$id)))
      stop("argument 'events_series' must be in the list of possible events_series.",call.=F)
      
      # transformer events_series en numerique
      if(is.character(events_series))
	events_series<-match(events_series,x@param$id)  
      
      # fonction de traçage
      .plot.cluster.one.1d(x,events_series=events_series,...)
      
    } else if(min(events_series)=="all") # toutes les séries d'évènement
    {
      
      # fonction de traçage
      .plot.cluster.all.1d(x,events_series=events_series)
      
    } else # plusieures les séries d'évènement (max 6)
    {   
      
      # vérifications sur events_series
      if(is.numeric(events_series) & max(events_series)>x@param$n_events_series)
        stop(paste("argument 'events_series' must be less than or equal to ",x@param$n_events_series,".",sep=""),call.=F)

      if(is.numeric(events_series) & min(events_series)<1)
        stop("argument 'events_series' must be greater than 0.",call.=F)      
     
      if(is.character(events_series) & any(is.na(match(events_series,x@param$id))))
        stop("argument 'events_series' must be in the list of possible events_series.",call.=F)
      
      # transformer events_series en numerique
      if(is.character(events_series))
        events_series<-match(events_series,x@param$id)
      
      # limiter la longeur de events_series
      if(length(events_series)>6)
        stop("argument 'events_series' must be of length less than 7.",call.=F)
      
      # fonction de traçage
      .plot.cluster.serie.1d(x,events_series=events_series)
    }
  }
  
  # --------------------------------------------
  # graphe 2d
  # --------------------------------------------  
  if(x@param$dimension=="2d")
  { 
      # fonction de traçage
      .plot.cluster.2d(x=x,map=map,indice=indice)
  }
  
  # --------------------------------------------
  # graphe 3d
  # --------------------------------------------  
  if(x@param$dimension=="3d")
  { 
      # fonction de traçage
      if(projection==FALSE)
	{ .plot.cluster.3d(x=x,indice=indice,sphere=sphere) } else
	{ .plot.cluster.3d.projection(x=x,indice=indice) }
  }
}

# --------------------------------
# méthode graphscan_plot
# --------------------------------
setGeneric(name="graphscan_plot",
	def=function(x,events_series=1,map=NULL,indice="cucala",sphere=TRUE,projection=FALSE,...){standardGeneric("graphscan_plot")}
)

setMethod(f="graphscan_plot",signature="graphscan",
  definition=function(x,events_series=1,map=NULL,indice="cucala",sphere=TRUE,projection=FALSE,...){ .plot.graphscan(x,events_series,map,indice,sphere,projection,...) })
  



