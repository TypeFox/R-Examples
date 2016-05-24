lsclass <-
function(type='numeric'){
	inlist<-ls(.GlobalEnv)
	classlist <- sapply(1:length(inlist),function(j) class(get(inlist[j])) )
	tnams<-sapply(1:length(inlist), function(j)type %in% classlist[[j]] )
	return(inlist[tnams])
}
