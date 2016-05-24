## Gain

setMethod("gain", signature=c(object="PSTf"), function(object, node, f, C) {
	node <- query(object, node, exact=TRUE, output="all")
	parent <- query(object, node.parent(node), output="all")

	seglist <- rownames(node@prob)

	res <- vector(mode="logical", length=length(seglist))
	names(res) <- seglist

	for (n in 1:length(seglist)) {
		if (f=="G2") {
			res[n] <- G2(node@prob[n,], parent@prob[seqglist[n],], C=C, N=node@n[n])
		} else if (f=="G1") {
			res[n] <- G1(node@prob[n,], parent@prob[seglist[n],], C=C)
		}
	}

	return(res)
}
)




