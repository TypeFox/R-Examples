.path_to_gb.dmp <-
function(package="ncbit"){
	path=paste(system.file(package=package), "data", sep="/")
#	pl=paths[[which(names(paths)=="fetch_genbank")]]
#	base_path=gsub("fetch_genbank.pl", "", pl)
	nd=paste(path, "nodes.dmp", sep="/")
	nm=paste(path, "names.dmp", sep="/")
	ncbi=paste(path, "ncbi.rda", sep="/")
	return(c(ncbi=ncbi, names=nm, nodes=nd, base=path))
}
