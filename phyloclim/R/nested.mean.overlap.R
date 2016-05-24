nested.mean.overlap <- function(phy, node, olap){
	
	# match ordering of phy and olap
	# ------------------------------
	id <- match(phy$tip.label, rownames(olap))
	olap <- olap[id, id]
	
	# get daughter nodes
	# ------------------
	d2 <- phy$edge[phy$edge[, 1] == node, 2]
	
	# get descendents of both daughter nodes
	# --------------------------------------
	C1 <- descendants(phy, d2[1])
	C2 <- descendants(phy, d2[2])
	
	# calculate mean overlap
	# ----------------------
	o <- 0
    for (j in C1){
        for (k in C2){
		    n <- nbConnectingNodes(phy, c(j, k))
	        o <- o + 0.5 ^ (n - 1) * olap[j, k]
	    }
    }
    o
}