# ---------------------------------------------------------------
# graphscan : version 1.1
# fonctions résumé des objets graphscan
#
# création   : 18/06/14
# version du : 18/06/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

# ---------------------------------------------------------------------------------
# fonction d'affichage des objets graphscan (utilisé par summary)
# ---------------------------------------------------------------------------------
.afficher_objet_graphscan_resume<-function(x)
{
  cat("Object of class graphscan - version 0.1\n")
  cat("\n")
  cat("Dimension                           : ",x@param$dimension,"\n")
  
  if(x@param$dimension=="1d")
  {
    if(!is.null(x@param$n_sequences)) 
      cat("Number of sequences in data         : ",x@param$n_sequences,"\n")
    cat("Number of events series used        : ",x@param$n_events_series,"\n") 
    if(!is.null(x@param$n)) 
      cat("Number of sites in sequences        : ",x@param$n,"\n")
    cat("Average of events number            : ",round(mean(x@param$n_events),1),"\n")
  }
  
  # affichage spécifique pour 2d et 3d
  if(x@param$dimension!="1d")
  {
    # calcul nombre de points avec et sans cas
    nb_point_cas<-table(x @data$x@data$cases>0)
    n<-nb_point_cas["TRUE"]
    n0<-nb_point_cas["FALSE"]
    # complexité nombre de distances cas-cas + nombre de distances cas-controles
    complexite<-(n^2-n)/2 + n*n0
    
    cat("Total number of cases               : ",x@param$n_events,"\n")
    cat(" Number of points with cases        : ",n,"\n")
    cat("Total number of controls            : ",x@param$n_controls,"\n") 
    cat(" Number of points without cases     : ",n0,"\n")
    cat("Number of distances to compute      : ",complexite,"\n")
  }  
  
  cat("\n")
  cat("Number of simulations               : ",x@param$n_simulation,"\n")
  cat("Significance threshold (alpha)      : ",x@param$alpha,"\n")
  cat("Type of analysis                    : ",x@param$cluster_analysis,"\n")  
  if(x@param$dimension=="1d") cat("User choice for the type of cluster : ",x@param$cluster_user_choice,"\n")
  cat("Concentration index                 : ",x@param$concentration_index,"\n\n")
 
}


setMethod(f="summary",signature="graphscan",
definition=function(object,...)
{
  # résumé pour objets 1D
  if(object@param$dimension=="1d")
  {
      # affichage des informations générales
      .afficher_objet_graphscan_resume(object)
      
      # extraction des informations sur les séries d'évènement
      res<-object@param$id
      res<-cbind(res,object@param$n_events)
      debut_fin<-t(sapply(object@param$normalisation_factor,FUN=function(x) x))
      res<-cbind(res,debut_fin)
      
      # mise au format data.frame
      res<-data.frame(res)
      names(res)<-c("events_series","n_events","start","end")
      
      # extraction si possible des résultats de l'analyse
      if(!is.null(object@cluster$cluster_1d_description)) 
      {
	cat("__________________________________________________________\n\n")
	cat("Analysis was performed on this object.\n")
	cat("__________________________________________________________\n")
	if(class(object@cluster$cluster_1d_description)=="matrix")
	{
	  res<-cbind(res,object@cluster$cluster_1d_description) 
	  names(res)[5:8]<-c("n_pos","n_neg","l_pos","l_neg")
	} else cat(object@cluster$cluster_1d_description,"\n\n")  
      }
      
     
  }
  
  
  # résumé pour objets 2D et 3D
  if(object@param$dimension=="2d" | object@param$dimension=="3d")
  {
      # affichage des informations générales
      .afficher_objet_graphscan_resume(object)
      
      if(length(object@cluster)>0)
      {
	  cat("__________________________________________________________\n\n")
	  cat("Analysis was performed on this object.\n")
	  cat("__________________________________________________________\n")
	  # cluster détecté avec les deux indices
	  if(!is.null(object@cluster$cluster_nd_cucala) & !is.null(object@cluster$cluster_nd_kulldorff))
	  {
	      res1<-object@cluster$cluster_nd_cucala@data[1,]
	      res1$name_index<-"cucala"
	      res2<-object@cluster$cluster_nd_kulldorff@data[1,]
	      res2$name_index<-"kulldorff"	      
	      res<-rbind(res1,res2)
	  }
	  
	  # cluster détecté avec indice de cucala seulement
	  if(!is.null(object@cluster$cluster_nd_cucala) & is.null(object@cluster$cluster_nd_kulldorff))
	  {
	      res1<-object@cluster$cluster_nd_cucala@data[1,]
	      res1$name_index<-"cucala"	 
	      res2<-object@cluster$cluster_nd_description[2]	      
	      res<-list(res1,res2)
	  }
	  
	  # cluster détecté avec indice de kulldorff seulement
	  if(is.null(object@cluster$cluster_nd_cucala) & !is.null(object@cluster$cluster_nd_kulldorff))
	  {
	      res1<-object@cluster$cluster_nd_description[1]	 
	      res2<-object@cluster$cluster_nd_kulldorff@data[1,]
	      res2$name_index<-"kulldorff"	      
	      res<-list(res1,res2)
	  }
	  
	  # aucun cluster détecté
	  if(is.null(object@cluster$cluster_nd_cucala) & is.null(object@cluster$cluster_nd_kulldorff))
	  {
	      res<-object@cluster$cluster_nd_description
	  }
      }
  }
  
  return(res)
})


# -------------------------------------------------------
# Autre nom pour la fonction 
# -------------------------------------------------------
# # associer la fonction générique S3 'summary' à la fonction 'summary.graphscan'
# setGeneric(name="summary.graphscan",def=function(object,...){standardGeneric("summary.graphscan")})
# 
# setMethod(f="summary.graphscan",signature="graphscan",
# definition=function(object,...){ summary(object=object,...) })




