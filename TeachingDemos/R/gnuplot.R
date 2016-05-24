gpenv <- new.env()
gpenv$gp <- numeric(0)
gpenv$gp.tempfiles <- character(0)

gp.open <- function(where='c:/progra~1/GnuPlot/bin/pgnuplot.exe'){
	gpenv$gp <<- pipe(where,'w')
	gpenv$gp.tempfiles <<- character(0)
	invisible(gpenv$gp)
}


gp.close <- function(pipe=gpenv$gp){
	cat("quit\n",file=pipe)
	close(pipe)
	if(exists('gpenv$gp.tempfiles')){
		unlink(gpenv$gp.tempfiles)
		gpenv$gp.tempfiles <- character(0)
	}
	gpenv$gp <<- numeric(0)
	invisible()
}

gp.send <- function(cmd='replot',pipe=gpenv$gp){
	cat(cmd, file=pipe)
	cat("\n",file=pipe)
	invisible()
}

gp.plot <- function(x,y,type='p',add=FALSE, title=deparse(substitute(y)),
		pipe=gpenv$gp){
	tmp <- tempfile()
	gpenv$gp.tempfiles <<- c(gpenv$gp.tempfiles, tmp)

	write.table( cbind(x,y), tmp, row.names=FALSE, col.names=FALSE )
	w <- ifelse(type=='p', 'points', 'lines')
	r <- ifelse(add, 'replot', 'plot')

	cat( paste(r," '",tmp,"' with ",w," title '",title,"'\n",sep=''),
		file=pipe)
	invisible()
}
