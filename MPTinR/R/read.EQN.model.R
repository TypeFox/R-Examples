
# this version of a function to read EQN files behaves like MultiTree, that is, it ignores the content of the first line and reads 
# all remaining lines!

.read.EQN.model <- function(model.filename) {
	parse.eqn <- function(x){
		branches <- unique(x[,2])
		l.tree <- length(branches)
		tree <- vector('expression', l.tree)
		for (branch in 1:l.tree) {
			tree[branch] <- parse(text = paste(x[x$V2 == branches[branch],"V3"], collapse = " + "))
		}
		tree
	}
	tmp.in <- read.table(model.filename, skip = 1, stringsAsFactors = FALSE)
	tmp.ordered <- tmp.in[order(tmp.in$V1),]
	tmp.spl <- split(tmp.ordered, factor(tmp.ordered$V1))
	tmp.spl <- lapply(tmp.spl, function(d.f) d.f[order(d.f[,2]),])
	model <- lapply(tmp.spl, parse.eqn)
	names(model) <- NULL
	model
}
