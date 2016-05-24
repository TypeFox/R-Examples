
# this version of a function to read EQN files behaves like HMMtree and GPT, that is, it only read the number of lines that are 
# indicated in the first line

.read.EQN.model.2 <- function(model.filename) {
	parse.eqn <- function(x){
		branches <- unique(x[,2])
		l.tree <- length(branches)
		tree <- vector('expression', l.tree)
		for (branch in 1:l.tree) {
			tree[branch] <- parse(text = paste(x[x$V2 == branches[branch],"V3"], collapse = " + "))
		}
		tree
	}
	n.lines <- read.table(model.filename, nrows = 1)[1,1]
	tmp.in <- read.table(model.filename, skip = 1, stringsAsFactors = FALSE, nrows = n.lines)
	tmp.ordered <- tmp.in[order(tmp.in$V1),]
	tmp.spl <- split(tmp.ordered, factor(tmp.ordered$V1))
	tmp.spl <- lapply(tmp.spl, function(d.f) d.f[order(d.f[,2]),])
	model <- lapply(tmp.spl, parse.eqn)
	names(model) <- NULL
	model
}
