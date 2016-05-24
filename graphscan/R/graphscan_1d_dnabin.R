# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction graphscan_1d : création des objets de classe graphscan 
# fonction pour data au format DNAbin
# création : 04/11/13
# version du : 04/11/13
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

setMethod(f="graphscan_1d",signature="DNAbin",
definition=function(data,format="fasta",events_series="all",id=NULL,n_simulation=199,cluster_analysis="both",normalisation_factor=NULL,alpha=0.05,cluster_user_choice="positive")
{
	# ----------------------------------
	# vérifications des paramètres
	# ----------------------------------
	.graphscan_1d_parametre_verification(format,events_series,id,n_simulation,
	                                    cluster_analysis,alpha,cluster_user_choice)
	
	
	# ----------------------------------
	# déterminer nombre de séquences et longeurs des séquences 
	# data doit être une matrice : avoir des dimensions
	# ----------------------------------
	if(is.null(dim(data))) stop("sequences must be of same length",call.=F)
	
	nb_sequence<-nrow(data)
	longeur<-ncol(data)
	
	# ----------------------------------
	# créer liste des séries d'évènement
	# et identifiants des séries d'évènement
	# ----------------------------------  
	liste_serie_evenement<-NULL
	
	# cas events_series = 'all'
	if(is.character(events_series)) 
	{
		nb_serie_evenement<-(nb_sequence**2-nb_sequence)/2
		liste_serie_evenement<-combn(x=1:nb_sequence,m=2)
	}
  

	 # cas events_series = list(A,B)
	if(is.list(events_series))
	{
	  if(length(events_series)!=2)
	      stop("'events_series' must be a list of length 2",call.=F)
	  message_erreur<-paste("indices in lists of 'events_series' must be comprise between 1 and",nb_sequence)
	  if(max(events_series[[1]])>nb_sequence | max(events_series[[2]])>nb_sequence)
	      stop(message_erreur,call.=F)
	  if(max(events_series[[1]])<1 | max(events_series[[2]])<1)
	      stop(message_erreur,call.=F)   
	  rm(message_erreur)
	  
	  # produit cartésien
	  liste_serie_evenement<-expand.grid(events_series[[1]],events_series[[2]])
	  
	  # éliminer comparaisons i vs i
	  iid<-liste_serie_evenement[,1]!=liste_serie_evenement[,2] 
	  liste_serie_evenement<-liste_serie_evenement[iid,]
	  
	  # éliminer doublon i vs j et j vs i
	  liste_serie_evenement<-t(apply(liste_serie_evenement,MARGIN=1,FUN=sort))
	  liste_serie_evenement<-unique(liste_serie_evenement)
	  liste_serie_evenement<-t(liste_serie_evenement)
      
	  # nombre de comparaisons
	  nb_serie_evenement<-dim(liste_serie_evenement)[2]
	}
	
	
	# identifiants des séquences
	if(is.null(id)) id<-"events_series_"
	creation_id_comparaison<-function(couple)
	{
	  x<-id_sequence[couple[1]]
	  y<-id_sequence[couple[2]]
	  res<-paste(id,x,"_vs_",y,sep="")
	  return(res)
	}
	id_sequence<-dimnames(data)[[1]]
	id_comparaison<-apply(X=liste_serie_evenement,MARGIN=2,FUN=creation_id_comparaison) 
  
	# ----------------------------------
	# créer liste des vecteurs de positions,
	# des vecteurs de positions normalisés et 
	# du vecteur du nombre d'évènements
	# le vecteur doit commencer par Zero
	# et se terminer par le facteur de
	# normalisation (longeur de la séquence +1)
	# ----------------------------------   
	determiner_position<-function(couple)
	{
	  x<-data[couple[1],]
	  y<-data[couple[2],]
	  position<-(1:length(x))[x!=y]
	  # ajour 0 et longeur+1 en début et fin du vecteur des positions
	  position<-c(0,position,as.integer(longeur+1))
	  
	  return(position)
	}
	
	calcul_nombre_evenement<-function(x)
	{
	      # -2 pour ne pas compter position 0 et longeur+1
	      return(length(x)-2)
	}
	
	position<-apply(X=liste_serie_evenement,MARGIN=2,FUN=determiner_position)
	if(!is.list(position)) position<-split(position,1)
	nb_evenement<-unlist(lapply(position,FUN=calcul_nombre_evenement))
	
	# ----------------------------------
	# bornes pour la normalisation
	# ----------------------------------
	normalisation<-rep(list(c(1,longeur+1)),times=nb_serie_evenement)
 
	# ----------------------------------
	# créer objet graphscan
	# ----------------------------------
	# liste des paramètres
	list_param<-list("id"=id_comparaison,"n_events_series"=as.integer(nb_serie_evenement),"dimension"="1d",
			"n"=as.integer(longeur),"normalisation_factor"=normalisation,
			"n_sequences"=as.integer(nb_sequence),"n_events"=as.integer(nb_evenement),
			"n_simulation"=as.integer(n_simulation),"cluster_analysis"=cluster_analysis,
			"alpha"=alpha,"cluster_user_choice"=cluster_user_choice,
			"concentration_index"="Cucala")                
			
			
	# liste des données
	list_data<-list("x"=position)
	
	res<-.create_graphscan(param=list_param,data=list_data)
	return(res)
})
