gnsc.restruc <-
function(pathway){
	path = unique(pathway)
	path.se.id = matrix(0,nrow=length(path),ncol=2)
	path.n=rep(0, length(path))
	for(k in 1:length(path)){
		tmp = (1:length(pathway))[pathway == path[k]]
		path.n[k] = length(tmp)
		path.se.id[k,1] = tmp[1]
		path.se.id[k,2] = tmp[path.n[k]]
	}
	output=list()
	output$path.n=path.n
	output$path.se.id=path.se.id
	return(output)
}
