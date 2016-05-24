Outetm <- function(data, v, ...) {
	# from state 1 to state 3
	id <- which(data[,v[1]] == data[,v[3]])
	idlen <- length(id)
	etm13 <- data.frame("id"=id, "from"=rep(1, times=idlen), "to"=rep(3, times=idlen), "entry"=rep(0, times=idlen), "exit"=data[id,v[1]], data[id,-v])
	cens <- which(data[,v[1]] == data[,v[3]] & data[,v[2]] == 0)
	for ( i in 1:length(cens) ) etm13$to[etm13$id == cens[i]] <- "cens"

	# from state 1 to state 2
	id <- which(data[,v[3]] > data[,v[1]])
	idlen <- length(id)
	etm12 <- data.frame("id"=id, "from"=rep(1, times=idlen), "to"=rep(2, times=idlen), "entry"=rep(0, times=idlen), "exit"=data[id,v[1]], data[id,-v])

	# from state 2 to state 3
	etm23 <- data.frame("id"=id, "from"=rep(2, times=idlen), "to"=rep(3, times=idlen), "entry"=data[id,v[1]], "exit"=data[id,v[3]], data[id,-v])
	cens <- which(data[,v[3]] > data[,v[1]] & data[,v[4]] == 0)
	for ( i in 1:length(cens) ) etm23$to[etm23$id == cens[i]] <- "cens"

	etmdata <- rbind(etm12, etm13, etm23)
	etmdata <- etmdata[order(etmdata$id),]
	row.names(etmdata) <- as.integer( 1:nrow(etmdata) )
	return(etmdata)
}
