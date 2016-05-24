tree.comp <-
function (phy1, phy2, method = "prop") {
	phy1 <- unroot(phy1)
	phy2 <- unroot(phy2)
	nphy1 <- length(phy1$tip.label)
	bp1 <- prop.part(phy1)
	bp1 <- lapply(bp1, function(xx) sort(phy1$tip.label[xx]))
	nphy2 <- length(phy2$tip.label)
	bp2.tmp <- prop.part(phy2)
	bp2 <- lapply(bp2.tmp, function(xx) sort(phy2$tip.label[xx]))
	bp2.comp <- lapply(bp2.tmp, function(xx) setdiff(1:nphy2, xx))
	bp2.comp <- lapply(bp2.comp, function(xx) sort(phy2$tip.label[xx]))
	q1 <- length(bp1)
	q2 <- length(bp2)
	p <- 0
        for (i in 1:q1) {
            for (j in 1:q2) {
                if (identical(bp1[[i]], bp2[[j]]) | identical(bp1[[i]], 
                  bp2.comp[[j]])) {
                  p <- p + 1
                  break
                }
            }
	}
	if(method == "shallow"){
	    shallow <- which(node.depth(phy1)[node.depth(phy1) > 1] <= median(node.depth(phy1)[node.depth(phy1) > 1]))
	    p2 <- 0
		for (i in shallow) {
		    for (j in 1:q2) {
			if (identical(bp1[[i]], bp2[[j]]) | identical(bp1[[i]], 
			  bp2.comp[[j]])) {
			  p2 <- p2 + 1
			  break
			}
		    }
		}
	dT <- p2/length(shallow)
	}
	if(method == "PH85") dT <- q1 + q2 - 2 * p
	if(method =="prop") dT <- p/q1
dT
}

