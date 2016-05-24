# ---------------------------------------------------------------
# graphscan : version 1.1
# fonctions cluster, .cluster_1d, .cluster_nd
# lancer les programmes C pour lancer les analyses de détection
# création : 23/10/13
# version du : 18/09/18
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

setGeneric(name="cluster",
	def=function(gr,n_simulation=NULL,cluster_analysis=NULL,memory_size=2000){standardGeneric("cluster")}
)


setMethod(f="cluster",signature="graphscan",
definition=function(gr,n_simulation=gr@param$n_simulation,cluster_analysis=gr@param$cluster_analysis,memory_size=2000)
{
	# -----------------------------------------------------
	# vérifications des paramètres de la fonction cluster
	# -----------------------------------------------------
	# argument n_simulation
	if(!is.numeric(n_simulation))
		stop("argument 'n_simulation' must be 'numeric'",call.=F)
	if(n_simulation<10 | n_simulation>10000)
		stop("argument 'n_simulation' must be comprise between 10 and 10000",call.=F)
	if(!is.numeric(memory_size) || memory_size<0)
		stop("argument 'memory_size' must be 'numeric' and >0",call.=F)
	
	# argument cluster_analysis
	if(is.na(match(cluster_analysis,c("both","positive","negative"))))
		stop("argument 'cluster_analysis' must be 'both', 'positive', 'negative'",call.=F)

	# modifier si besoin n_simulation et cluster_analysis de l'objet gr
	if(n_simulation!=gr@param$n_simulation) gr@param$n_simulation<-n_simulation
	if(cluster_analysis!=gr@param$cluster_analysis) gr@param$cluster_analysis<-cluster_analysis
	
	# lancer l'analyse
	if(gr@param$dimension == "1d")
	{
		gr<-.cluster_1d(gr,n_simulation,cluster_analysis)
	}
	else
	{
		gr<-.cluster_nd(gr,n_simulation,cluster_analysis,memory_size)
	}
	return(gr)
})


# ----------------------------------------------------
# fonction pour détection 1d
# ----------------------------------------------------

.cluster_1d<-function(gr,n_simulation=gr@param$n_simulation,cluster_analysis=gr@param$cluster_analysis)
{
	resultat<-NULL
	cat("in progress ...\n")
	for(i in 1:gr@param$n_events_series)
	{
	    # -----------------------------------------------------
	    # récupérer les paramètres de l'objet gr 
	    # -----------------------------------------------------
	    # nombre d'évènements
	    nb_evenement<-gr@param$n_events[i]
	    
	    # bornes pour la normalisation
	    normalisation_debut<-gr@param$normalisation_factor[[i]][1]
	    normalisation_fin<-gr@param$normalisation_factor[[i]][2]
	    
	    # vecteurs des positions des évènements et indicatrice de normalisation
	    vecteur_evenement<-as.numeric(gr@data$x[[i]])
	    
	    # seuil de la significativité des clusters
	    alpha<-gr@param$alpha
	    
	    # paramètre de serie_evenement de taille des agrégats
	    theta<-10
	    
	    # nombre de simulations pour le calcul de la significativité
	    # géré en fonction du paramètre n_simulation de la fonction cluster
	    if(n_simulation==gr@param$n_simulation) nb_simulation<-gr@param$n_simulation else nb_simulation<-n_simulation
	    
	    # choix du type de détection
	    # géré en fonction du paramètre cluster_analysis de la fonction cluster
	    if(cluster_analysis==gr@param$cluster_analysis) choix_detection_r_string<-gr@param$cluster_analysis else choix_detection_r_string<-cluster_analysis
	    choix_detection<-match(choix_detection_r_string,c("positive","negative","both")) # codage numérique 1, 2, 3.
	    
	    # dans le cas de détection des clusters positifs et negatif en même temps cluster_analysis='both'
	    # dans le cas de deux agregats, un positif et un négatif, autant significatif l'un que l'autre
	    # choix du type de cluster
	    choix_type_agregat<-match(gr@param$cluster_user_choice,c("negative","positive","random")) # codage numérique 1, 2, 3.
	    
	    # --------------------------------------------------------------------------------
	    # exécuter la recherche de cluster par la fonction C  'detection_multiple_dagregat'
	    # pour données 1d
	    # --------------------------------------------------------------------------------
	    # la fonction C 'detection_multiple_dagregat' renvoie les résultats sous la forme :
	    # start,end,index,pvalue,positivity,id_cluster,id_serie
	    
	    res<-NULL
	    cat("events series : ",i,"/",gr@param$n_events_series,"\n")
	    res<-.Call("detection_multiple_dagregat",nb_evenement,normalisation_debut,
				    normalisation_fin,vecteur_evenement,alpha,theta,
				    nb_simulation,choix_detection,choix_type_agregat)
	    
            # gérer résultat vide
            if(ncol(res)==0) res<-matrix(NA,nrow=6,ncol=1)
	    
	    res<-rbind(res,rep(i,times=dim(res)[2]))
	    resultat<-rbind(resultat,t(res))  
	}
	
	# --------------------------------------------------------------------------------
	# traitements du tableau de résultats
	# --------------------------------------------------------------------------------
	if(!is.na(resultat[1][1])) # vérifier présence d'un résultat
	  {
	    res<-.cluster_1d_traitement_tableau(resultat,gr)
	    resultat_brut<-res[[1]]
	    resultat_traite<-res[[2]]
	    description_serie_evenement<-res[[3]]
	  } else
	  {
	    resultat_brut<-NULL
	    resultat_traite<-NULL
	    # message si aucun cluster n'est détecté
	    description_serie_evenement<-paste("No cluster detected at level alpha = ",gr@param$alpha,sep="")
	  }
	
	# ajout résultats de l'analyse à l'objet gr
	gr@cluster[["cluster_1d"]]<-resultat_traite
	gr@cluster[["cluster_1d_raw"]]<-resultat_brut
	gr@cluster[["cluster_1d_description"]]<-description_serie_evenement
	
	# renvoyer l'objet gr
	return(gr)
}


# ----------------------------------------------------
# fonction pour détection nd
# ----------------------------------------------------

.cluster_nd<-function(gr,n_simulation=gr@param$n_simulation,cluster_analysis=gr@param$cluster_analysis,memory_size=2000)
{
	resultat<-NULL

	# nombre de simulations pour le calcul de la significativité
	# en fonction du paramètre n_simulation de la fonction cluster
	if(n_simulation==gr@param$n_simulation)
	  nb_simulation<-gr@param$n_simulation else nb_simulation<-n_simulation

	# dimension des coordonnées (format numérique)
	dimension<-match(gr@param$dimension,c("1d","2d","3d"))
		
	# coordonnées des points au format x1,y1,x2,y2,...
	coordonnees_r<-as.vector(t(gr@data$x@coords))
	
	# ensemble des identifiants des points
	id_r<-1:nrow(gr@data$x)
	id_r<-as.integer(id_r)
	
	# nombre de cas et controles par point
	cas_r<-as.numeric(gr@data$x@data$cases)
	controle_r<-gr@data$x@data$controls
	
	# nombre total de points (cas + controles)
	# nb_point<-gr@param$n_events+gr@param$n_controls
	# modification David
	nb_point<-nrow(gr@data$x)

	# recherche des clusters
	resultat = .Call("detection_cluster",nb_point,dimension,n_simulation,id_r,coordonnees_r,controle_r,cas_r,as.integer(memory_size))

	# mise en forme des résultats au format SpatialPointsDataFrame
	# mise à l'exponentielle de l'indice de kulldorff
	nb_point_cucala<-resultat[3,1]
	if(nb_point_cucala>0)
		cucala_data<-data.frame(id=resultat[7:(6+nb_point_cucala),1], index=rep(resultat[1,1],times =nb_point_cucala),radius=rep(resultat[2,1],times=nb_point_cucala),pvalue=rep(resultat[4,1],times=nb_point_cucala),n_controls=rep(resultat[5,1],times=nb_point_cucala),n_cases=rep(resultat[6,1],times=nb_point_cucala))

	nb_point_kulldorf<-resultat[3,(dimension+2)]
	if(nb_point_kulldorf>0)
		kulldorff_data<-data.frame(id=resultat[7:(6+nb_point_kulldorf),(dimension+2)], index=rep(resultat[1,(dimension+2)],times=nb_point_kulldorf),radius=rep(resultat[2,(dimension+2)],times=nb_point_kulldorf),pvalue=rep(resultat[4,(dimension+2)],times=nb_point_kulldorf),n_controls=rep(resultat[5,(dimension+2)],times=nb_point_kulldorf),n_cases=rep(resultat[6,(dimension+2)],times=nb_point_kulldorf))
	
	
	# description des résultats
	if(nb_point_cucala>0 && nrow(cucala_data)>0 && cucala_data[1,4]<=gr@param$alpha)
	{
		cucala<-SpatialPointsDataFrame(coords=resultat[7:(6+nb_point_cucala),2:(dimension+1)],data= cucala_data)
	  	description_cucala<-paste("Cucala index: ",sprintf(fmt="%10.3e",cucala_data[1,2])," - pvalue: ",round(cucala_data[1,4],3),sep="")
	  	description_cucala<-paste(description_cucala," \nn_cases: ",cucala_data[1,6]," - n_controls: ",cucala_data[1,5],sep="")
	  	description_cucala<-paste(description_cucala," - radius: ",round(cucala_data[1,3],2),sep="")
	}else{
	  	description_cucala<-paste("No cluster detected with Cucala index at level alpha = ",gr@param$alpha,sep="")
	  	cucala<-NULL
	}
	
	# description des résultats
	if(nb_point_kulldorf>0 && nrow(kulldorff_data)>0 && kulldorff_data[1,4]<=gr@param$alpha)
	{
		kulldorff<-SpatialPointsDataFrame(coords=resultat[7:(6+nb_point_kulldorf),(dimension+3):(2*dimension+2)],data= kulldorff_data)
	  	description_kulldorff<-paste("Kulldorff index: ",sprintf(fmt="%10.3e",kulldorff_data[1,2])," - pvalue: ",round(kulldorff_data[1,4],3),sep="")
	  	description_kulldorff<-paste(description_kulldorff," \nn_cases: ",kulldorff_data[1,6]," - n_controls: ",kulldorff_data[1,5],sep="")
	  	description_kulldorff<-paste(description_kulldorff," - radius: ",round(kulldorff_data[1,3],2),sep="")
	}else{
	  	description_kulldorff<-paste("No cluster detected with Kulldorff index at level alpha = ",gr@param$alpha,sep="")
	  	kulldorff<-NULL
	}
	
	# ajout des résultats à l'objet gr
	gr@cluster<-list()
	gr@cluster[["cluster_nd_cucala"]]<-cucala
	gr@cluster[["cluster_nd_kulldorff"]]<-kulldorff
	gr@cluster[["cluster_nd_description"]]<-c(description_cucala,description_kulldorff)
	return(gr)
}
