# ---------------------------------------------------------------
# graphscan : version 1.1
# fonction graphscan_get_number_proc
# récupérer le nombre de processeurs sur la machine
# création : 16/01/2015
# Unité Epidémiologie Animale (UR346)
# Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
# ---------------------------------------------------------------

.graphscan_get_number_proc <- function(){
	tmp<-.C('get_number_proc',res=as.integer(0))
	return(tmp$res)
}