terminal.sisters <- function(phy){

	x <- sapply(1:length(phy$tip.label), sister, phy = phy)
	len <- sapply(x, length)
	x[len > 1] <- NA
	x <- unlist(x)
	x <- cbind(1:length(phy$tip.label), x)
	x <- x[!is.na(x[, 2]), ]
	x <- cbind(phy$tip.label[x[, 1]], phy$tip.label[x[, 2]])
	unique(t(apply(x, 1, sort)))
}