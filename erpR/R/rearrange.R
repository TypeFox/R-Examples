rearrange <-
function(deps, oth=NULL, dataset, name.dep="Dep", name.newvar="New_Var") {
	dep=as.vector(as.matrix(dataset[,deps]))
	others=NULL
	if (length(oth)>0){
		for (i in 1:dim(dataset[,deps])[2]){
			others=rbind(others, dataset[,oth])
		}
	}
	nuova.var=NULL
	for (i in 1:length(names(dataset[,deps]))) {
		nuova.var=c(nuova.var,rep(names(dataset[,deps])[i], dim(dataset[,deps])[1]))
		}
	dat=as.data.frame(cbind(dep,others))
	names(dat)[1]=paste(name.dep)
	dat[,paste(name.newvar)]=as.factor(nuova.var)
	return(dat)
}
