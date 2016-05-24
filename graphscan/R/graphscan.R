# ---------------------------------------------------------------
# graphscan : version 1.1
# définition de la classe graphscan
# fonction 'print' pour afficher les informations d'un objet de classe graphscan
#
# création : 23/10/13
# version du : 05/05/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

# listes pour prototypes de la classe graphscan
# (valeurs par défaut)
list_param<-list("id"=character(),"n_events_series"=integer(1),"dimension"=character(1),"n_sequences"=integer(),
                 "n"=integer(),"normalisation_factor"=NULL,"n_events"=integer(),
                 "n_controls"=integer(),"n_simulation"=as.integer(199),"cluster_analysis"="both",
                 "alpha"=0.05,"cluster_user_choice"="positive","concentration_index"="both")
                 
list_data<-list("x"=list(),"point_events"=NULL,"point_controls"=NULL)

list_cluster<-list()

# fonction de création d'objets de classe graphscan
.create_graphscan<-setClass(Class="graphscan",
slots = c(param="list",data="list",cluster="list"),
prototype = list(param=list_param,data=list_data,cluster=list_cluster))

# effacer listes pour prototypes
rm(list_param,list_data,list_cluster)

# ---------------------------------------------------------------------------------
# fonction d'affichage des objets graphscan (utilisé par print et show)
# ---------------------------------------------------------------------------------
.afficher_objet_graphscan<-function(x)
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
  
  # affichage spécifique à la 2d et 3d
  if(x@param$dimension!="1d")
  {
    # calcul nombre de points avec et sans cas
    nb_point_cas<-table(x@data$x@data$cases>0)
    n<-nb_point_cas["TRUE"]
    n0<-nb_point_cas["FALSE"]
    
    cat("Total number of cases               : ",x@param$n_events,"\n")
    cat(" Number of points with cases        : ",n,"\n")
    cat("Total number of controls            : ",x@param$n_controls,"\n") 
    cat(" Number of points without cases     : ",n0,"\n")
  }  
  
  cat("\n")
  cat("Number of simulations               : ",x@param$n_simulation,"\n")
  cat("Significance threshold (alpha)      : ",x@param$alpha,"\n")
  cat("Type of analysis                    : ",x@param$cluster_analysis,"\n")  
  if(x@param$dimension=="1d") cat("User choice for the type of cluster : ",x@param$cluster_user_choice,"\n")
  cat("Concentration index                 : ",x@param$concentration_index,"\n")
  
  # affichage si possible des résultats de l'analyse
  if(x@param$dimension=="1d" & !is.null(x@cluster$cluster_1d_description)) 
  {
  cat("__________________________________________________________\n\n")
  cat("Analysis output for the five first events series.\n")
  cat("__________________________________________________________\n")
  if(class(x@cluster$cluster_1d_description)=="matrix") print(head(x@cluster$cluster_1d_description,n=5)) else cat(x@cluster$cluster_1d_description,"\n\n")  
  }
  
  if(x@param$dimension!="1d" & length(x@cluster)>0)
  {
  cat("__________________________________________________________\n\n")
  cat("Analysis output for Cucala and Kulldorff index.\n")
  cat("__________________________________________________________\n")
  cat(x@cluster$cluster_nd_description[1],"\n\n")
  cat(x@cluster$cluster_nd_description[2],"\n\n")
  }
}

# ---------------------------------------------------------------------------------
# fonction 'print.graphscan' pour afficher les informations d'un objet de classe graphscan
# ---------------------------------------------------------------------------------
# associer la fonction générique S3 'print' à la fonction .afficher_objet_graphscan
# pour un objet de classe 'graphscan'
setMethod(f="print",signature="graphscan",definition=.afficher_objet_graphscan)

# ---------------------------------------------------------------------------------
# fonction 'show.graphscan' pour afficher les informations d'un objet de classe graphscan
# ---------------------------------------------------------------------------------
setMethod(f="show",signature="graphscan",definition=function(object){.afficher_objet_graphscan(object)})
