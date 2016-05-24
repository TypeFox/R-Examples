# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction graphscan_1d : création des objets de classe graphscan 
# fonction pour data contenant une liste de vecteurs de positions d'évènement
# création : 06/11/13
# version du : 06/11/13
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

 
setMethod(f="graphscan_1d",signature="list",
definition=function(data,id=NULL,n_simulation=199,cluster_analysis="both",
                    normalisation_factor,alpha=0.05,cluster_user_choice="positive")
{
  # ----------------------------------
  # vérifications des paramètres
  # ----------------------------------
  events_series<-"all"
  format<-"fasta"
  .graphscan_1d_parametre_verification(format,events_series,id,n_simulation,cluster_analysis,
                                      alpha,cluster_user_choice)
                                       
  # ----------------------------------
  # identifiants des vecteurs de positions
  # on utilise les noms des élèments de la liste
  # ----------------------------------  
  # nombre de vecteurs de positions
  if(length(data)>0) 
    nb_vecteur<-length(data) else stop("'data' list is empty")
  # identifiants
  if(is.null(id)) id<-"vector_"
  if(is.null(names(data)))
    id_events_series<-paste(id,1:nb_vecteur,sep="") else id_events_series<-names(data)

  # ----------------------------------
  # vérification absence des données manquantes (NA)
  # ----------------------------------
  verification_donnee_manquante<-function(x) return(any(is.na(x)))
  
  test_donnee_manquante<-any(unlist(lapply(data,FUN=verification_donnee_manquante)))
  if(test_donnee_manquante) 
    stop("any vector contain NA.\n",call.=F)
  
  # ----------------------------------
  # nombre d'évènement par vecteur
  # ----------------------------------
  nb_evenement<-unlist(lapply(data,FUN=length))
  test_nb_evenement<-any(nb_evenement==0)
  if(test_nb_evenement)
    stop("one vector is empty.",call.=F)
 
  # ----------------------------------
  # tri des vecteurs
  # ---------------------------------- 
  data<-lapply(data,FUN=sort)
 
  # ----------------------------------
  # normalisation des données
  # ----------------------------------
  # vérification de la normalisation
  if(is.null(normalisation_factor))
      stop("'normalisation_factor' must be defined",call.=F)
  if(!is.list(normalisation_factor))
    stop("'normalisation_factor' must be a list of vector.",call.=F)
  if(length(normalisation_factor)!=nb_vecteur)
      stop("list 'normalisation_factor' must be of same length than vector of positions",call.=F)    

  
  # vérifier que les bornes sont ordonnées
  test_borne<-any(unlist(lapply(normalisation_factor,FUN=function(x) return(x[2]>x[1]))))
  if(!test_borne)
    stop("in 'normalisation_factor' the bounds must be ordered.",call.=F)
   
  
  # vérifier que les positions sont bien comprises entre les bornes
  for(i in 1:nb_vecteur)
  {
    res<-data[[i]]
    a<-normalisation_factor[[i]][1]
    b<-normalisation_factor[[i]][2]
    
    # vérification de la normalisation
    if(min(res)<a | max(res)>b)
      stop(paste("normalisation failed for vector ",i,sep=""),call.=F)
    
    # ajout de la borne inférieure au début 
    # et supérieure à la fin du vecteur si non comprise
    if(min(data[[i]])>a) data[[i]]<-c(a,data[[i]])
    if(max(data[[i]])<b) data[[i]]<-c(data[[i]],b)
    
  }
  
  # ----------------------------------
  # vérification des ex-aequo
  # ----------------------------------
  verification_exaequo<-function(x) return(any(duplicated(x)))
  #elimination_exaequo<-function(x) return(unique(x))
  
  test_exaequo<-any(unlist(lapply(data,FUN=verification_exaequo)))
  if(test_exaequo) 
    stop("any vector contain duplicated positions.\n",call.=F)
  
 
  # ----------------------------------
  # créer objet graphscan
  # ----------------------------------
  # liste des paramètres
  list_param<-list("id"=id_events_series,"n_events_series"=as.integer(nb_vecteur),"dimension"="1d",
                   "n"=NULL,"n_sequences"=NULL,
                   "n_events"=nb_evenement,
                   "normalisation_factor"=normalisation_factor,
                   "n_simulation"=as.integer(n_simulation),"cluster_analysis"=cluster_analysis,
                   "alpha"=alpha,"cluster_user_choice"=cluster_user_choice,
                   "concentration_index"="Cucala")
                   
  # liste des données
  list_data<-list("x"=data)
  
  res<-.create_graphscan(param=list_param,data=list_data)
  return(res)
  
})