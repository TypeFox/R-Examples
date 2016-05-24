# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction cluster_1d_array_treatment : traitement du tableau de résultat 1d 
# renvoie : un tableau de données brutes, un tableau de données traitées (prise ne compte des
# des clusters inclus dans les autres et donc du découpage en segements) et d'un tableau de
# description des clusters
#
# création : 14/04/14
# version du : 14/04/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------
	
.cluster_1d_traitement_tableau<-function(resultat,gr)
{

	  # recoder identifiants des clusters détectés à partir de 1
	  resultat[,6]<-resultat[,6]+1
	  
	  # création des tableau resultat_brut, resultat_traite et description_serie_evenement
	  description_serie_evenement<-NULL
	  resultat_brut<-NULL	
	  resultat_traite<-NULL
	  
	  
	  # liste des séries d'évènement
	  liste_serie_evenement<-unique(resultat[,7])
	  nb_serie_evenement<-length(liste_serie_evenement)
	  
	  # ----------------------------------------------------
	  # calculs des coordonnées des segments représentant les clusters
	  # et indice de cucala, en tenant compte des clusters inclus
	  # ----------------------------------------------------
	  # début et fin des segements
	  xleft_total<-NULL
	  xright_total<-NULL
	  
	  # indice, pvaleur et positivité des segments
	  indice_total<-NULL
	  positivite_total<-NULL 
	  p_valeur_total<-NULL
	  
	  # nombre de segments nécessaires pour chaque cluster
	  # (prise en compte des clusters inclus)
	  n_segment_total<-NULL 
	
	  # identifiants des segments (leur identifiant de cluster) et de la série
	  id_segment_total<-NULL
	  id_serie_total<-NULL
	  
	  ########################################################### section parallélisée avec snowfall
	  # nombre maximum de processeurs
	  nb_proc_max<-.graphscan_get_number_proc()
	  
	  # nombre de processeur à utiliser
	  nb_proc<-min(nb_serie_evenement,nb_proc_max)
	  
	  # fonction extraction des coordonnées segments
	  # à utiliser pour list.graphscan_1d_parameter_verificatione_k une séries d'évènement
	  
	  calcul_segment<-function(liste_k){
	  	  # initialisation des vecteurs de sortie
		  n_segment_total<-NULL
		  xleft_total<-NULL
		  xright_total<-NULL
		  indice_total<-NULL
		  id_segment_total<-NULL
		  positivite_total<-NULL
		  id_serie_total<-NULL
		  p_valeur_total<-NULL
			
		  for(k in liste_k) # pour chaque élèment de la liste
		  {
			# extraire les données pour la série d'évènement en cours
			d<-resultat[resultat[,7]==k,]
			if(is.null(dim(d))) d<-matrix(d,nrow=1) # garder format matrix
			
			# nombre de clusters pour la séries d'évènement en cours
			nb_cluster<-nrow(d)
			if(is.na(d[1,1])) nb_cluster<-0
			
			# -------------------------------------
			# détecter cluster inclus
			# -------------------------------------			
			inclusion<-list(NULL)
			
			if(nb_cluster>1) # si on plus d'un cluster
			{
			      # trier les clusters par ordre de détection
			      id<-sort(d[,6],index.return=T)
			      d<-d[id$ix,]
			      if(is.null(dim(d))) d<-matrix(d,nrow=1) # garder format matrix
			      
			      # détecter cluster inclus dans un autre 
			      # pour créer segments et calculer longeurs
			      # si un cluster A contient un autre cluster B
			      # on cherche l'identifiant du cluster inclus B	
			      for(j in 1:nb_cluster)
			      {
				test_inclusion<-d[j,1]<=d[-j,1] & d[j,2]>=d[-j,2]
				inclusion[[d[j,6]]]<-d[-j,6][test_inclusion]
			      }

			      # dans la liste des clusters inclus dans le cluster A 
			      # enlever les clusters déjà inclus dans un autre cluster que A
			      for(j in 2:nb_cluster) for(l in 1:(j-1))
			      {
				id<-match(inclusion[[l]],inclusion[[j]])
				id<-id[!is.na(id)]	  
				if(length(id)>0) inclusion[[j]]<-inclusion[[j]][-id]  
			      }
			      

			      # -------------------------------------
			      # tri sur les identifiants des clusters
			      # -------------------------------------
			      # traiter les clusters dans l'ordre de détection
			      id<-sort(d[,6],index.return=T)
			      d<-d[id$ix,]
			      if(is.null(dim(d))) d<-matrix(d,nrow=1) # garder format matrix
			      
			      # -------------------------------------
			      # calculs des coordonnées des segments
			      # dessinant les clusters
			      # -------------------------------------
			      # ces segments doivent prendre en compte les
			      # clusters inclus dans les autres et ainsi découper les segments
			      
			      # traitement par cluster
			      for(j in 1:nb_cluster)
			      { 
				# coordonnées et identifiants des segments pour le cluster en cours
				xleft<-xright<-NULL
				indice<-NULL   
				id_segments<-NULL
				positivite<-NULL
				id_serie<-NULL
				p_valeur<-NULL
				
				# cas sans clusters inclus
				if(length(inclusion[[j]])==0) 
				{
				  
				  n_segment<-1
				  xleft<-d[j,1]
				  xright<-d[j,2]
				  indice<-d[j,3]
				  id_segment<-j
				  positivite<-d[j,5]
				  id_serie<-d[j,7]
				  p_valeur<-d[j,4]
				  
				} else # présence de clusters inclus
				{
				  
				  n_segment<-length(inclusion[[j]])+1
				  indice<-rep(d[j,3],times=n_segment)
				  id_segment<-rep(j,times=n_segment)
				  positivite<-rep(d[j,5],times=n_segment)
				  id_serie<-rep(d[j,7],times=n_segment)
				  p_valeur<-rep(d[j,4],times=n_segment)
				  n_segment<-rep(n_segment,times=n_segment)
				  
				  # ouverture du 1er segment pour le cluster en cours
				  xleft<-d[j,1]
				  
				  # extraire clusters inclus dans le cluster en cours
				  dd<-d[inclusion[[j]],]
				  if(is.null(dim(dd))) dd<-matrix(dd,nrow=1) # garder format matrix		      
				    
				  # tri sur les débuts des clusters inclus
				  id<-sort(dd[,1],index.return=T)
				  dd<-dd[id$ix,]
				  if(is.null(dim(dd))) dd<-matrix(dd,nrow=1) # garder format matrix
				  
				  # mise en forme vectorielle
				  dd<-as.vector(t(dd[,1:2]))
				  # extraire une coordonnée sur deux à partir de 2
				  id<-seq(from=2,to=length(dd),by=2)
				  xleft<-c(d[j,1],dd[id])
				  xright<-dd[id-1]
				  
				  # fermeture du dernier segment pour le cluster en cours
				  xright<-c(xright,d[j,2])
				}     
			    
				# remplir vecteurs généraux par 
				# les coordonnées et le nombre de segments
				n_segment_total<-c(n_segment_total,n_segment)
				xleft_total<-c(xleft_total,xleft)
				xright_total<-c(xright_total,xright)
				indice_total<-c(indice_total,indice) 
				id_segment_total<-c(id_segment_total,id_segment)
				positivite_total<-c(positivite_total,positivite)
				id_serie_total<-c(id_serie_total,id_serie)
				p_valeur_total<-c(p_valeur_total,p_valeur)
			      }
		      # fin condition nb_cluster>1
		      } else if(nb_cluster==1) # cas nb_cluster==1
		      {
		      
				n_segment_total<-c(n_segment_total,1)
				xleft_total<-c(xleft_total,d[1,1])
				xright_total<-c(xright_total,d[1,2])
				indice_total<-c(indice_total,d[1,3]) 
				id_segment_total<-c(id_segment_total,1)
				positivite_total<-c(positivite_total,d[1,5])
				id_serie_total<-c(id_serie_total,d[1,7])
				p_valeur_total<-c(p_valeur_total,d[1,4])
				
		      } else # cas sans cluster
		      {
		      
				n_segment_total<-c(n_segment_total,NA)
				xleft_total<-c(xleft_total,NA)
				xright_total<-c(xright_total,NA)
				indice_total<-c(indice_total,NA) 
				id_segment_total<-c(id_segment_total,NA)
				positivite_total<-c(positivite_total,NA)
				id_serie_total<-c(id_serie_total,k)
				p_valeur_total<-c(p_valeur_total,NA)
				
		      }
		      
		} # fin boucle sur longeurlist_k
		
		return(list(n_segment_total,xleft_total,xright_total,indice_total,
		            id_segment_total,id_serie_total,positivite_total,p_valeur_total))
	  } # fin fonction
	  

	  # distribution des charges entre processeurs 
	  list_k<-list(1)
	  if(nb_serie_evenement>nb_proc)
	  {
	      charge<-as.numeric(cut(1:nb_serie_evenement,breaks=nb_proc))
	      list_k<-split(1:nb_serie_evenement,f=charge)
	  } else if(nb_serie_evenement>1)
	  {
	      list_k<-as.list(1:nb_serie_evenement)
	  }
	  
	  # exécution en parallele ou non
	  if(length(list_k)>1)
	  {
	      # initialisation pour nb_proc processeurs
	      sfInit(parallel=TRUE,cpus=nb_proc,socketHosts=NULL)
	      # exporter la matrice resultat aux processeurs
	      sfExport("resultat")
	      # executions de la fonction calcul_segment
	      res<-sfLapply(x=list_k,fun=calcul_segment)
	      # arreter parallélisation
	      sfStop() 
	   } else
	   {
	      res<-calcul_segment(1)
	   }
	  
	  # assemblage des résultats. Le format de res dépend du nombre
	  # du séries d'évènement présentes.
	  if(nb_serie_evenement>1)
	  {
	      # format list de list
	      n_segment<-unlist(lapply(res,FUN=function(x) x[[1]]))
	      xleft<-unlist(lapply(res,FUN=function(x) x[[2]]))
	      xright<-unlist(lapply(res,FUN=function(x) x[[3]]))
	      index<-unlist(lapply(res,FUN=function(x) x[[4]]))
	      id_segment<-unlist(lapply(res,FUN=function(x) x[[5]]))
	      id_serie<-unlist(lapply(res,FUN=function(x) x[[6]]))
	      positivity<-unlist(lapply(res,FUN=function(x) x[[7]]))
	      pvalue<-unlist(lapply(res,FUN=function(x) x[[8]]))
	  } else
	  {
	      # format list
	      n_segment<-res[[1]]
	      xleft<-res[[2]]
	      xright<-res[[3]]
	      index<-res[[4]]
	      id_segment<-res[[5]]
	      id_serie<-res[[6]]
	      positivity<-res[[7]]
	      pvalue<-res[[8]]
	  }
	  
	  
	  ############################################################# fin section parallélisée
	  
	  # ------------------------------------------------------
	  # création du tableau des résultats traités
	  # ------------------------------------------------------
	  resultat_traite<-cbind(xleft,xright,index,pvalue,positivity,n_segment,id_segment,id_serie)
	  
	  # ajouter la variable longeur
	  resultat_traite<-cbind("length"=resultat_traite,resultat_traite[,2]-resultat_traite[,1])
	  
	  # éliminer segments de longeur nulle pour séries avec cluster
	  id<-resultat_traite[,9]!=0
	  id[is.na(id)]<-TRUE 
	  resultat_traite<-resultat_traite[id,]
          if(is.null(dim(resultat_traite))) resultat_traite<-matrix(resultat_traite,nrow=1) # garder format matrix		      
	  
	  # éliminer doublons
	  resultat_traite<-unique(resultat_traite)
	  if(is.null(dim(resultat_traite))) resultat_traite<-matrix(resultat_traite,nrow=1) # garder format matrix
		
	  # noms des lignes et des colonnes
	  rownames(resultat_traite)<-1:nrow(resultat_traite)
	  colnames(resultat_traite)<-c("xleft","xright","index","pvalue","positivity","n_segment","id_segment","id_serie","length")
		
	  # ------------------------------------------------------
	  # création du tableau descriptif des clusters
	  # ------------------------------------------------------			     	     
	  # facteur de normalisation (séries sur le même segment)
	  normalisation<-unlist(lapply(gr@param$normalisation_factor,FUN=function(x) x[2]-x[1]))
	  		     
	  # calcul des proportions des longeurs	
	  # des clusters positifs et négatifs série par série	  
	  longeur<-tapply(resultat_traite[,9],INDEX=list(resultat_traite[,8],resultat_traite[,5]),FUN=sum)
	  longeur[is.na(longeur)]<-0
	  
	  # ajout de la colonne manquante (si tous les clusters sont positifs ou négatifs)
	  if(ncol(longeur)==1) if(colnames(longeur)=="0") longeur<-cbind(longeur,0) else longeur<-cbind(0,longeur)
	  
	  # normalisation
	  longeur[,1]<-longeur[,1]/normalisation*100
	  longeur[,2]<-longeur[,2]/normalisation*100 
	  longeur<-round(longeur,2)
	  
	  # nombre de clusters positifs et négatifs
	  nb_cluster<-tapply(resultat_traite[,9],INDEX=list(resultat_traite[,8],resultat_traite[,5]),FUN=length)	
	   # ajout de la colonne manquante (si tous les clusters sont positifs ou négatifs)
	  if(ncol(nb_cluster)==1) if(colnames(nb_cluster)=="0") nb_cluster<-cbind(nb_cluster,0) else nb_cluster<-cbind(0,nb_cluster)
	  
	  # mettre cluster positif en premier
	  # on ajoutte [,1:2] (à gauche) pour garder le format matrix dans
	  # le cas d'une seule série d'évènement à analyser.
	  longeur[,1:2]<-longeur[,2:1]
	  nb_cluster[,1:2]<-nb_cluster[,2:1]
	  
	  # création tableau et identifiant de la serie_evenement dans le nom de la ligne
	  description_serie_evenement<-cbind(nb_cluster,longeur)
	  description_serie_evenement[is.na(description_serie_evenement)]<-0
	  colnames(description_serie_evenement)<-c("n_pos","n_neg","l_pos","l_neg")
	  rownames(description_serie_evenement)<-gr@param$id		      

	  
	  # ------------------------------------------------------
	  # modification du tableau des résultats bruts
	  # ------------------------------------------------------	  
	  resultat_brut<-resultat
	  colnames(resultat_brut)<-c("xleft","xright","index","pvalue","positivity","id_segment","id_serie")
	  
	  # renvoyer résultats
	  return(list(resultat_brut,resultat_traite,description_serie_evenement))
}