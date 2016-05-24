Predict.Plot <- function(model, pred.var, ...,
		type='response', add=FALSE, plot.args=list(),
		n.points=100, ref.val,
		ref.col='green', ref.lty=1, data ) {

	x2 <- list(...)

	if(missing(pred.var)) pred.var <- names(x2)[1]

	if(is.character(plot.args)) plot.args <- eval(parse(text=plot.args))

	getdata <- function(model) {
		if ('data' %in% names(model)) return(model$data)
		tmpcall <- model$call
		tmpcall[[1]] <- as.name('glm')
		model <- eval(tmpcall)
		model$data
	}


	if( pred.var %in% names(x2) ) {
		if (length(x2[[pred.var]]) > 1) {
			tmp.x <- seq( min(x2[[pred.var]]), max(x2[[pred.var]]),
				length.out=n.points)
		} else {
			if( missing(data) ) data <- getdata(model)
			ref.val <- x2[[pred.var]]
			tmp.x <- seq( min(data[[pred.var]]), max(data[[pred.var]]),
				length.out=n.points)
		}
	} else {
		if( missing(data) ) data <- getdata(model)
		tmp.x <- seq( min(data[[pred.var]]),  max(data[[pred.var]]),
			length.out=n.points)
	}

	x2[[pred.var]] <- tmp.x
	x <- as.data.frame(x2)

	yhat <- predict(model, x, type=type)

	if(add){
		plot.args$x <- x[[pred.var]]
		plot.args$y <- yhat
		do.call(lines, plot.args)
	} else {
		nms <- names(plot.args)
		plot.args$x=x[[pred.var]]
		plot.args$y=yhat
		if( !( 'ylab' %in% nms ) ) plot.args$ylab='Predicted Value'
		if( !( 'xlab' %in% nms ) ) plot.args$xlab=pred.var
		if( !( 'type' %in% nms ) ) plot.args$type='l'
		do.call(plot, plot.args)
	}


	if(!missing(ref.val)){
		tmp.x <- list(...)
		tmp.x[[pred.var]] <- ref.val
		yhat <- predict(model, as.data.frame(tmp.x), type=type)
		usr <- par('usr')
		lines( c(ref.val,ref.val,usr[1]), c(usr[3],yhat,yhat),
			col=ref.col, lty=ref.lty)
	}
}




TkPredict <- function(model, data, pred.var, ...){
	if( missing(data) ){
		if( class(model)[1] == 'lm' ){
			tmpcall <- model$call
			tmpcall[[1]] <- as.name('glm')
			model2 <- eval(tmpcall)
		} else {
			model2 <- model
		}
		data <- model2$data
	}

	tr <- delete.response( terms(model) )
	x <- get_all_vars(tr, data)

	if(missing(pred.var)) pred.var <- names(x)[1]

	lst <- list()

	lst$pred.var <- list('radiobuttons',values=names(x), init=pred.var)
        lst[[2]] <- list()
	for ( v in names(x) ) {
		tmp.x <- x[[v]]
		if( is.factor(tmp.x) ) {
			lvls <- levels(tmp.x)
			if(length(lvls) < 11 ) {
				lst[[2]][[v]] <- list('radiobuttons', values=lvls,
					init=lvls[1] )
			} else {
				lst[[2]][[v]] <- list('Entry', init=lvls[1])
			}
		} else {
			tmp.min <- min(tmp.x)
			tmp.max <- max(tmp.x)
			tmp.med <- median(tmp.x)
			lst[[2]][[v]] <- list('slider',from=tmp.min, to=tmp.max,
				init=tmp.med, resolution=signif( (tmp.max-tmp.min)/100, 2 ) )
		}
	}

        lst[[3]] <- list()
	lst[[3]]$plot.args <- list( 'entry', init='list()' )
	lst[[3]]$type <- list('entry', init='response')

	cl <- as.call( substitute( Predict.Plot(model) ) )

	eval(substitute(tkexamp( cl, lst, plotloc='left' )))

}




