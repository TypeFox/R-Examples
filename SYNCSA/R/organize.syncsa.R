organize.syncsa <-
function (comm, traits, dist.spp, envir) 
	{
	if (missing(comm)=="TRUE"){
		stop("\n Community not fount\n")
		}
	if (is.null(colnames(comm))){
		stop("\n Erro in column names of comm\n")
		}
	if (is.null(rownames(comm))){
		stop("\n Erro in row names of comm\n")
		}	
    if (!missing(traits) == "TRUE") {
    	if (is.null(colnames(traits))){
		stop("\n Erro in column names of traits\n")
		}
		if (is.null(rownames(traits))){
		stop("\n Erro in row names of traits\n")
		}
		match.names<-match(colnames(comm), rownames(traits))
		if(sum(is.na(match.names))>0)
		{
		stop("\n There are species from community data that are not on traits matrix\n")
		}
		traits<-as.matrix(traits[match.names,])
		}
	if (!missing(dist.spp) == "TRUE") {
		if (is.null(colnames(dist.spp))){
		stop("\n Erro in column names of dist.spp\n")
		}
		if (is.null(rownames(dist.spp))){
		stop("\n Erro in row names of dist.spp\n")
		}
		match.names<-match(colnames(comm), colnames(dist.spp))
		if(sum(is.na(match.names))>0)
		{
		stop("\n There are species from community data that are not on phylogenetic distance matrix\n")
		}
	dist.spp<-as.matrix(dist.spp[match.names,match.names])
	}
	if (!missing(envir) == "TRUE") {
		if (is.null(colnames(envir))){
		stop("\n Erro in column names of envir\n")
		}		
		if (is.null(rownames(envir))){
		stop("\n Erro in row names of envir\n")
		}
		match.names<-match(rownames(comm), rownames(envir))
		if(sum(is.na(match.names))>0)
		{
		stop("\n The are community that are not on environmental matrix\n")
		}
		envir<-as.matrix(envir[match.names,])   
    }
    if (missing(traits)=="TRUE"){
	traits<-NULL
		}
	if (missing(dist.spp)=="TRUE"){
	dist.spp<-NULL
		}
	if (missing(envir)=="TRUE"){
	envir<-NULL
		}		
    return(list(community=as.matrix(comm),traits=traits,dist.spp=dist.spp,environmental=envir))
}

