gp.splot <- function(x, y, z, add= FALSE, title=deparse(substitute(z)),
	pipe=gpenv$gp, datafile=tempfile()) {

	tmp <- datafile
	gpenv$gp.tempfiles <- c(gpenv$gp.tempfiles, tmp)
	tmp2 <- data.frame(x=x, y=y, z=z)
	tmp2 <- tmp2[ order(x,y), ]
	tmp3 <- split(tmp2, tmp2$x)
	con <- file(tmp, open='w')
	sapply( tmp3, function(d) {
		write.table( d, con, row.names=FALSE, col.names=FALSE )
		cat( "\n", file=con )
		} )
	close(con)
	cat( ifelse(add, "replot", "splot"), " '", tmp, "' title '",
		title, "'\n", sep="", file=pipe )
	invisible()
}

