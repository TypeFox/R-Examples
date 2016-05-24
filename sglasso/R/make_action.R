make_action <- function(theta){
	th_nm <- rownames(theta)
	np <- dim(theta)[2]
	id <- abs(theta) > 0
	mdl.list <- apply(id, 2, function(z) th_nm[z])
	if(np > 1) {
		action <- rep("", length = np)
		action[1] <- paste(mdl.list[[1]], collapse = " + ")
		for(j in 2:np) action[j] <- setDiff(mdl.list[[j]],mdl.list[[j-1]])
	} else action <- paste(drop(mdl.list), collapse = " + ")
	action
}
