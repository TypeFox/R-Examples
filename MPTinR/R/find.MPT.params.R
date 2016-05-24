
# .find.MPT.params <- function(model) {
	# inobjects <- function(ex) {
		# split.exp1 <- strsplit(as.character(ex),"[[:space:]]")
		# split.exp2 <- sapply(split.exp1, strsplit, split = "[()+*-]")
		# return(sort(unique(grep("[[:alpha:]]",unlist(split.exp2), value = TRUE))))
		# }
	# tmp <- sapply(model,inobjects)
	# return(unique(sort(unique(unlist(tmp)))))
# }

.find.MPT.params <- function(model) {
	tmp <- lapply(model,all.vars)
	return(unique(sort(unique(unlist(tmp)))))
}

