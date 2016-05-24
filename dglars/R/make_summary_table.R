make_summary_table <- function(x,k,complexity){
	n <- dim(x$X)[1]
	np <- x$np
	action <- x$action
	b <- x$beta
	dev <- x$dev
	g <- x$g
	if(complexity=="df")	compl <- x$df
	else compl <- gdf(x)
	gof <- dev + k * compl
	rank.gof <- rank(gof)
	best <- rank.gof==1
	b.gof <- b[,best]
	b.gof <- b.gof[abs(b.gof)>0]
	mark <- rep("  ",np)
	mark[best] <- "<-"
	rank.gof <- paste(rank.gof,mark)
	tbl <- data.frame(Sequence=action,g=g,Dev=dev,Complexity=compl,gof=gof,Rank=rank.gof)
	list(table=tbl,b.gof=b.gof)
}
