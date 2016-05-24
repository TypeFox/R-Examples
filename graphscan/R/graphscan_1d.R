# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction graphscan_1d : création des objets de classe graphscan 
# fonction générique et vérifications des paramètres
#
# création : 23/10/13
# version du : 20/03/14
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------
setGeneric(name="graphscan_1d",
    def=function(data,format="fasta",events_series="all",id=NULL,n_simulation=199,cluster_analysis="both",normalisation_factor=NULL,alpha=0.05,cluster_user_choice="positive"){standardGeneric("graphscan_1d")}
)


# --------------------------------------------
# fonction de vérifications des paramètres 1d
# --------------------------------------------
.graphscan_1d_parametre_verification<-function(format,events_series,id,n_simulation,cluster_analysis,alpha,cluster_user_choice)
{
  # argument format
  if(is.na(match(format,c("interleaved","sequential","clustal","fasta"))))
    stop("argument 'format' must be 'interleaved', 'sequential', 'clustal' or 'fasta'",call.=F)
  
  # argument events_series
  if(!(is.character(events_series) | is.list(events_series)))
    stop("argument 'events_series' must be 'all' or 'list'",call.=F)
  if(is.character(events_series)) if(events_series!="all")
    stop("argument 'events_series' must be 'all' or 'list'",call.=F)
  if(is.list(events_series) & length(events_series)!=2)
    stop("argument 'events_series' must be 'list' of length 2",call.=F)
    
  # argument id
  if(!(is.null(id) | is.character(id)))
    stop("argument 'id' must be 'NULL' or 'character'",call.=F)
  if(is.character(id) & length(id)!=1)
    stop("argument 'id' must be 'character' of length 1",call.=F)
  
  # argument n_simulation
  if(!is.numeric(n_simulation))
    stop("argument 'n_simulation' must be 'numeric'",call.=F)
  if(n_simulation<10 | n_simulation>10000)
    stop("argument 'n_simulation' must be comprise between 10 and 10000",call.=F)
  
  # argument cluster_analysis
  if(is.na(match(cluster_analysis,c("both","positive","negative"))))
    stop("argument 'cluster_analysis' must be 'both', 'positive' or 'negative'",call.=F)
    
  # argument alpha seuil de la significativité des clusters
  if(!is.numeric(n_simulation))
    stop("argument 'alpha' must be 'numeric'",call.=F)
  if(alpha<0 | alpha>1)
    stop("argument 'alpha' must be comprise between 0 and 1",call.=F)
  
  # dans le cas de détection des clusters positifs et negatif en même temps cluster_analysis='both'
  # dans le cas de deux agregats, un positif et un négatif, autant significatif l'un que l'autre
  # choix du type de cluster
  # argument cluster_user_choice
  if(is.na(match(cluster_user_choice,c("negative","positive","random"))))
    stop("argument 'cluster_user_choice' must be 'positive', 'negative' or 'random'",call.=F)
  
  
  return(NULL)
}
