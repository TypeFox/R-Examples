plot2script <- function(file='clipboard'){
	con <- file(file)
	open(con, open='a')
	tmp <- recordPlot()[[1]]
	for (i in seq(along.with=tmp)){
		fn <- tmp[[i]][[1]]
		args <- tmp[[i]][[2]]
		fns <- deparse(fn)
		m <- sub('^.*"(.*)".*$', '\\1', fns, perl=TRUE)
		c2 <- as.list(c(m,args))
		tmp2 <- do.call('call',c2)
                tmp3 <- match.call(get(m), call=tmp2)
                if(tmp3[[1]] == 'box'){
                    tmp3$which <- c("plot", "figure", "inner",
                                    "outer")[ tmp3$which ]
                }
		dput(tmp3, file=con)
	}
	close(con)
}

oldzoomplot <- function( xlim, ylim=NULL ){
	xy <- xy.coords(xlim,ylim)
	xlim <- range(xy$x)
	ylim <- range(xy$y)

	tmp <- recordPlot()[[1]]
	for(i in seq(along=tmp)){
		fn <- tmp[[i]][[1]]
		alst <- as.list(tmp[[i]][[2]])
		tmp2 <- all.equal( '.Primitive("locator")', deparse(fn) )
		if(is.logical(tmp2) && tmp2){
			next
		}
		tmp2 <- all.equal( '.Primitive("plot.window")', deparse(fn) )
		if(is.logical(tmp2) && tmp2) {
			alst[[1]] <- xlim
			alst[[2]] <- ylim
		}

		do.call(fn, alst)
	}
}



zoomplot <- function( xlim, ylim=NULL ){
	xy <- xy.coords(xlim,ylim)
	xlim <- range(xy$x)
	ylim <- range(xy$y)

	tmp <- recordPlot()
	for(i in seq(along=tmp[[1]])){
		fn <- tmp[[1]][[i]][[2]][[1]]
		if(fn$name == 'C_plot_window') {
			tmp[[1]][[i]][[2]][[2]] <- xlim
			tmp[[1]][[i]][[2]][[3]] <- ylim
		}
	}
        replayPlot(tmp)
}


