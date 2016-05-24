"parOKgroups" <-
function(clu,
parOKaddParam #list of additional parameters, at lest k and groups
){
	isTRUE(all(cut(clu,c(0,cumsum(parOKaddParam$k)),labels =FALSE)==parOKaddParam$groups))
	
}
