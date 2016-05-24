seqetm<-function(seq, method="transition", use.labels=TRUE, sep=">", bp="", ep="end"){

	statl <- alphabet(seq)#seqstatl(seq)
	nr <- attr(seq, "nr")
	has.nr <- any(seq==nr)
	if (has.nr) {
		statl <- c(statl, nr)
	}
	nbstat <- length(statl)
	tevent <- matrix(nrow=nbstat, ncol=nbstat)
	rownames(tevent) <- statl
	colnames(tevent) <- statl
	alphabet <- statl
	if (use.labels && inherits(seq, "stslist")) {
		#label<-alphabet(seq)
		label <- attr(seq, "labels")
		if (has.nr) {
			label <- c(label, nr)
		}
		if(length(label)==length(alphabet)){
			alphabet <- label
		}
		else if(length(label)>0){
			warning("Length of the labels and of the alphabet are not equal")
		}
	}
	if(any(grepl(",", alphabet))){
		warning(" [!] Alphabet and/or state labels should not contain commas ',' which are reserved for separating multiple events of a same transition!\n")
	}
	for(i in 1:nbstat){
		for(j in 1:nbstat){
			if(i==j){
				tevent[i,j] <- alphabet[[i]]
			}else{
				if(method=="transition"){
					tevent[i,j] <- paste(alphabet[[i]], alphabet[[j]], sep=sep)
				}else if(method == "period"){
					tevent[i,j] <- paste(ep, alphabet[[i]], ",", bp, alphabet[[j]], sep="")
				}else if(method == "state"){
					tevent[i,j] <- alphabet[[j]]
				}
			}
		}
	}
	return(tevent)

}