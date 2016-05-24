make_action <- function(beta,np){
	vn <- rownames(beta)
	id <- abs(beta)>0
	mdl.list <- apply(id,2,function(z) vn[z])
	if(np==1) action <- "Int."
	else {
		action <- rep("",length=np)
		for(j in 2:np) action[j-1] <- setDiff(mdl.list[[j]],mdl.list[[j-1]])
	}
	action
}
