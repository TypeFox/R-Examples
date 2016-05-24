`getListIdentity` <-
function(mat,n.cat){
	if(any(is.na(mat)))
		tmp.na<-!is.na(mat)
	else
		tmp.na<-TRUE
	listIdentity<-vector("list",n.cat)
	for(i in 1:n.cat)
		listIdentity[[i]]<-tmp.na & mat==i
	listIdentity
}

