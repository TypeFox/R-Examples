# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction graphscan_nd : création des objets de classe graphscan 
# en dimension 2d, 3d ou 4d 
# création : 16/12/13
# version du : 06/05/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

setGeneric(name="graphscan_nd",
	def=function(data,field_cases=NULL,field_controls=NULL,n_simulation=199,alpha=0.05){standardGeneric("graphscan_nd")}
)

# --------------------------------------------
# fonction de vérifications des paramètres nd
# --------------------------------------------

.graphscan_nd_parametre_verification<-function(data,n_simulation,alpha,cluster_analysis)
{
	# argument data
	if(class(data)!="SpatialPointsDataFrame")
		stop("argument 'data' must be of class 'SpatialPointsDataFrame'",call.=F)    
	
	# argument n_simulation
	if(!is.numeric(n_simulation))
		stop("argument 'n_simulation' must be 'numeric'",call.=F)
	if(n_simulation<10 | n_simulation>10000)
		stop("argument 'n_simulation' must be comprise between 10 and 10000",call.=F)
	
	# argument cluster_analysis
	if(is.na(match(cluster_analysis,c("both","positive","negative"))))
		stop("argument 'cluster_analysis' must be 'both', 'positive' or 'negative'")
	
	# argument alpha seuil de la significativité des clusters
	if(!is.numeric(alpha))
		stop("argument 'alpha' must be 'numeric'",call.=F)
	if(alpha<0 | alpha>1)
		stop("argument 'alpha' must be comprise between 0 and 1",call.=F)
	
	return(NULL)
}


# --------------------------------------------
# fonction graphscan_nd
# --------------------------------------------
setMethod(f="graphscan_nd",signature="SpatialPoints",
definition=function(data,field_cases=NULL,field_controls=NULL,n_simulation=199,alpha=0.05)
{
	# ----------------------------------
	# vérifications des paramètres
	# ----------------------------------
	cluster_analysis<-"positive" # seule méthode implémentée
	.graphscan_nd_parametre_verification(data,n_simulation,alpha,cluster_analysis)
	
	# ----------------------------------
	# dimension
	# ----------------------------------
	dimension<-paste(ncol(data@coords),"d",sep="")
	
	# ----------------------------------
	# récupération ou création des nombres de cas par point
	# ----------------------------------
	# si aucun nom de champ pour le nombres de cas par point
	# on cherche si un champ 'cases' existe.
	if(is.null(field_cases))
	{
		if(match("cases",names(data),nomatch=0)==0)
		{
		  stop("no field found for the number of cases in data",call.=F)
		} else # le champ 'cases' existe
		{	
		  case<-data@data$cases
		}
	} else # le nom du champ cases est précisé dans field_cases
	{	
		id_cases<-match(field_cases,names(data),nomatch=0)
		# vérifier l'existence du champ
		if(id_cases==0)
		  stop(paste("field for the number of cases '",field_cases,"' not exist in data",sep=""),call.=F)
		case<-data@data[,id_cases]
		
		# créer champ cases
		data@data$cases<-case
	}
	
	# ----------------------------------
	# récupération ou création des nombres de controle par point
	# ----------------------------------
	# si aucun nom de champ pour le nombres de controle par point
	# on cherche si un champ 'controls' existe.
	if(is.null(field_controls))
	{
		if(match("controls",names(data),nomatch = 0) == 0)
		{
		  stop("no field found for the number of controls in data",call.=F)
		} else # le champ'controls' existe
		{
		    controle<-data@data$controls
		}
	} else # le nom du champ cases est précisé dans field_controls
	{
		id_control<-match(field_controls,names(data),nomatch=0)
		# vérifier l'existence du champ
		if(id_control==0)
		  stop(paste("field for the number of controls '",field_controls,"' not exist in data",sep=""),call.=F)
		controle<-data@data[,id_control]
		
		# créer champ controls
		data@data$controls<-controle
	}
	
	
	# ----------------------------------
	# récupérer nombre de cas 
	# et nombre de controles et vérifications
	# ----------------------------------
	nb_evenement<-sum(case)
	nb_controle<-sum(controle)
	
 	if(nb_evenement<5)
 		stop("the number of cases must be equal or greater than 5",call.=F)
 	if(nb_controle<50)
 		stop("the number of controls must be equal or greater than 50",call.=F)
	
	
	# ----------------------------------
	# créer objet graphscan nd
	# ----------------------------------
	# liste des paramètres
	list_param<-list("dimension"=dimension,
			 "n_events"=as.integer(nb_evenement),"n_controls"=as.integer(nb_controle),
			 "n_simulation"=as.integer(n_simulation),"cluster_analysis"=cluster_analysis,
			 "alpha"=alpha,"concentration_index"="Cucala and Kulldorff")
	
	# liste des données
	list_data<-list("x"=data)
	
	res<-.create_graphscan(param=list_param,data=list_data)
	return(res)
})

