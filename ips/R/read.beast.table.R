read.beast.table <- function(file, digits = 2){
	phy <- read.beast(file, digits = digits)
	int <- phy$Nnode
	stats <- phy[-(1:4)]
	M <- matrix(unlist(stats), nrow = int, byrow = FALSE)
	colnames(M)  <- names(stats)
	tips <- length(phy$tip.label)
	node <- (tips + 1):(tips + int)
	M <- cbind(node, M)
	M
}