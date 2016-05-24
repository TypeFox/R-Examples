pop.trajectories.plotAll <- function(pop.pred, 
									output.dir=file.path(getwd(), 'pop.trajectories'),
									output.type="png", expression=NULL, verbose=FALSE,  ...) {
	# plots pop trajectories for all countries
	if(!is.null(expression) && !grepl('XXX', expression, fixed=TRUE))
		stop('Expression must contain a mask "XXX" to be used as country code.')
	bayesTFR:::.do.plot.all.country.loop(pop.pred$countries[,'code'], meta=NULL, country.table=pop.pred$countries,
					output.dir=output.dir, func=.trajectories.plot.with.mask.replacement, output.type=output.type, 
					file.prefix='pop.plot', plot.type='population graph', verbose=verbose, pop.pred=pop.pred, 
					expression=expression, ...)
}

.trajectories.plot.with.mask.replacement <- function(pop.pred, country, expression=NULL, plotfun='pop.trajectories.plot', ...) {
	expr.arg <- NULL
	country.arg <- NULL
	if(!is.null(expression))
		expr.arg <- gsub('XXX', as.character(country), expression, fixed=TRUE)
	else country.arg <- country
	do.call(plotfun, list(pop.pred, country=country.arg, expression=expr.arg, ...))
}

pop.trajectories.plot <- function(pop.pred, country=NULL, expression=NULL, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  sum.over.ages=FALSE,
								  half.child.variant=FALSE,
								  nr.traj=NULL, typical.trajectory=FALSE, main=NULL,
								  dev.ncol=5, lwd=c(2,2,2,2,1), col=c('black', 'red', 'red', 'blue', '#00000020'),
								  show.legend=TRUE, ann=par('ann'), ...
								  ) {
	# lwd is a vector of 5 line widths for: 
	#	1. observed data, 2. median, 3. quantiles, 4. half child variant, 5. trajectories
	if (is.null(country) && is.null(expression)) {
		stop('Argument "country" or "expression" must be given.')
	}
	if(length(lwd) < 5) {
		llwd <- length(lwd)
		lwd <- rep(lwd, 5)
		lwd[(llwd+1):5] <- c(2,2,2,2,1)[(llwd+1):5]
	}
	if(length(col) < 5) {
		lcol <- length(col)
		col <- rep(col, 5)
		col[(lcol+1):5] <- c('black', 'red', 'red', 'blue', '#00000020')[(lcol+1):5]
	}

	if(!is.null(country)) {
		country <- get.country.object(country, country.table=pop.pred$countries)
		if(is.null(country$code)) stop('Country not available.')
	}
	if(sum.over.ages || age[1]=='psr' || !is.null(expression))
		do.pop.trajectories.plot(pop.pred, country, expression=expression, pi=pi, sex=sex, age=age,
									half.child.variant=half.child.variant, nr.traj=nr.traj,
									typical.trajectory=typical.trajectory,
									main=main, lwd=lwd, col=col,
									show.legend=show.legend, ann=ann, ...)
	else {
		all.ages <- pop.pred$ages
		if(age[1]=='all') age <- 1:20
		age.labels <- get.age.labels(pop.pred$ages)
		if(is.null(main)) {
			main <- country$name
			sex <- match.arg(sex) 
			if(sex != 'both') main <- paste(main, ': ', sex, sep='')
		}
		age.labels <- get.age.labels(pop.pred$ages)
		cur.mgp <- par('mgp')
		cur.oma <- par('oma')
		cur.mar <- par('mar')
		nplots <- length(age)
		if (nplots < dev.ncol) {
        	ncols <- nplots
			nrows <- 1
        } else {
			ncols <- dev.ncol
			nrows <- ceiling(nplots/dev.ncol)
        }		
		par(mfrow=c(nrows,ncols),  oma = c(0, 0, 2, 0))
		par(mar=c(2,2,1,0.4)+0.1, mgp=c(1,0.3,0))
		for(iage in age) {
			do.pop.trajectories.plot(pop.pred, country, pi=pi, sex=sex, age=iage,
									half.child.variant=half.child.variant, nr.traj=nr.traj,
									typical.trajectory=typical.trajectory,
									xlab='', ylab='', main=age.labels[iage], cex.main=0.9, 
									lwd=lwd, col=col, show.legend=show.legend, ann=ann, ...)
		}
		if(ann) mtext(main, line = 0.5, outer = TRUE)
		par(mgp=cur.mgp, mar=cur.mar, oma=cur.oma)
	}
}

do.pop.trajectories.plot <- function(pop.pred, country=NULL, expression=NULL, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  half.child.variant=FALSE,
								  nr.traj=NULL, typical.trajectory=FALSE,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='', ylab='Population projection', main=NULL, 
								  lwd=c(2,2,2,2,1), col=c('black', 'red', 'red', 'blue', '#00000020'),
								  show.legend=TRUE, ann=par('ann'), add=FALSE, adjust=FALSE, adj.to.file=NULL, ...
								  ) {

	sex <- match.arg(sex)
	if(!is.null(adj.to.file)) {
		adjust <- FALSE
		if (is.null(country))
			stop('Argument "country" must be given for the adjustment.')
	}
	reload.traj.if.needed <- !adjust
	if(!is.null(expression)) {
		trajectories <- get.pop.trajectories.from.expression(expression, pop.pred, nr.traj, 
										typical.trajectory=typical.trajectory, adjust=adjust, adj.to.file=adj.to.file, adj.country=country$code)
		if(missing(xlim) || (!missing(xlim) && min(xlim) < min(pop.pred$proj.years-4)))
			pop.observed.all <- get.pop.observed.from.expression(expression, pop.pred)
		else {
			pop.observed.all <- NA
			names(pop.observed.all) <- min(pop.pred$proj.years)-5
		}
		reload.traj.if.needed <- FALSE
	} else {
		trajectories <- get.pop.trajectories(pop.pred, country$code, sex, age, nr.traj, 
										typical.trajectory=typical.trajectory, adjust=adjust)
		pop.observed.all <- get.pop.observed(pop.pred, country$code, sex=sex, age=age)
	}
	#stop('')
	cqp <- list()
	for (i in 1:length(pi))
		cqp[[i]] <- get.pop.traj.quantiles(trajectories$quantiles, pop.pred, country$index, country$code, 
										trajectories=trajectories$trajectories, pi=pi[i], sex=sex, age=age, reload=reload.traj.if.needed)
	
	obs.not.na <- !is.na(pop.observed.all)
	pop.observed.idx <- if(sum(obs.not.na)==0) length(pop.observed.all) else which(obs.not.na)
	x1 <- as.integer(names(pop.observed.all)[pop.observed.idx])
	x2 <- if(!is.null(dimnames(trajectories$trajectories)) && !is.null(dimnames(trajectories$trajectories)[[1]])) 
				as.numeric(dimnames(trajectories$trajectories)[[1]])
			else as.numeric(dimnames(pop.pred$quantiles)[[3]])
	y1 <- pop.observed.all[pop.observed.idx]
	if(is.null(xlim)) xlim <- c(min(x1, x2), max(x1, x2))
	if(is.null(ylim)) 
		ylim <- c(min(y1, if (!is.null(trajectories$trajectories))
							trajectories$trajectories[,trajectories$index]
						  else NULL, 
						  sapply(cqp, min, na.rm=TRUE), na.rm=TRUE), 
				  max(y1, if (!is.null(trajectories$trajectories))
				  			trajectories$trajectories[,trajectories$index] else NULL, 
				  		  sapply(cqp, max, na.rm=TRUE), na.rm=TRUE))
	if(is.null(main)) {
		if(!is.null(expression)) main <- expression
		else {
			main <- country$name 
			if(sex != 'both') main <- paste(main, ': ', sex, sep='')
			if(age[1] == 'psr') main <- paste(main, ' (Potential Support Ratio)', sep='')
			else {
				if(age[1] != 'all') {
					age.labels <- get.age.labels(pop.pred$ages[age], collapsed=TRUE)
					main <- paste(main, ' (Age ', paste(age.labels, collapse=','), ')', sep='')
				}
			}
		}
	}
	# plot historical data: observed
	if(!add) {
		if(missing(ylab) && !is.null(expression)) ylab <- ''
		plot(x1, y1, type=type, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
			panel.first = grid(), lwd=lwd[1], col=col[1], ann=ann, ...)
	}
	else lines(x1, y1, type=type, lwd=lwd[1], col=col[1])
	# plot trajectories
	if(!is.null(trajectories$index) && !is.null(trajectories$trajectories)) {
		for (i in 1:length(trajectories$index)) {
			lines(x2, trajectories$trajectories[,trajectories$index[i]], type='l', col=col[5], lwd=lwd[5])
		}
	}
	# plot median
	pop.median <- get.pop.traj.quantiles(trajectories$quantiles, pop.pred, country$index, country$code, 
										trajectories=trajectories$trajectories, q=0.5, sex=sex, age=age, reload=reload.traj.if.needed)
	lines(x2, pop.median, type='l', col=col[2], lwd=lwd[2]) 
	
	# plot given CIs
	lty <- 2:(length(pi)+1)
	for (i in 1:length(pi)) {		
		if (!is.null(cqp[[i]])) {
			lines(x2, cqp[[i]][1,], type='l', col=col[3], lty=lty[i], lwd=lwd[3])
			lines(x2, cqp[[i]][2,], type='l', col=col[3], lty=lty[i], lwd=lwd[3])
		}
	}
	legend <- c('median', paste(pi, '% PI', sep=''))
	lty <- c(1, lty)
	lwds <- c(lwd[2], rep(lwd[3], length(pi)))
	cols <- c(col[2], rep(col[3], length(pi)))
	if(sum(obs.not.na)>0) {
		legend <- c(legend, 'observed')
		lty <- c(lty, 1)
		lwds <- c(lwds, lwd[1])
		cols <- c(cols, col[1])
	}
	if (half.child.variant && !is.null(trajectories$half.child)) {
		lty <- c(lty, max(lty)+1)
		llty <- length(lty)
		lines(x2, trajectories$half.child[,1], type='l', col=col[4], lty=lty[llty], lwd=lwd[4])
		lines(x2, trajectories$half.child[,2], type='l', col=col[4], lty=lty[llty], lwd=lwd[4])
		legend <- c(legend, '+/- 0.5 child')
		cols <- c(cols, col[4])
		lwds <- c(lwds, lwd[4])
	}
	if(show.legend && ann)
		legend('topleft', legend=legend, lty=lty, bty='n', col=cols, lwd=lwds)
}

pop.trajectories.table <- function(pop.pred, country=NULL, expression=NULL, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  half.child.variant=FALSE, ...) {
	do.pop.trajectories.table(pop.pred, country=country, expression=expression, pi=pi,
								  sex=sex, age=age, half.child.variant=half.child.variant, ...)					  	
}

do.pop.trajectories.table <- function(pop.pred, country=NULL, expression=NULL, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  half.child.variant=FALSE, adjust=FALSE, adj.to.file=NULL, ...) {
	if (is.null(country)  && is.null(expression)) 
		stop('Argument "country" or "expression" must be given.')
				
	if(!is.null(adj.to.file)) {
		adjust <- FALSE
		if (is.null(country))
			stop('Argument "country" must be given for the adjustment.')
	}
	
	if(!is.null(country))
		country <- get.country.object(country, country.table=pop.pred$countries)
	
	max.age.idx <- length(pop.pred$ages)
	sex <- match.arg(sex)
	quant <- NULL
	if (age[1]=='all') age.idx <- 1:max.age.idx
	else {
		if(all(is.element(1:max.age.idx, age))) age.idx <- 1:max.age.idx
		else age.idx <- unique(age)
	}
	lage <- length(age.idx)
	if(!is.null(expression)) {
		trajectories <- get.pop.trajectories.from.expression(expression, pop.pred, adjust=adjust, adj.to.file=adj.to.file, adj.country=country$code, ...)
		pop.observed <- get.pop.observed.from.expression(expression, pop.pred)
		quant <- NULL
	} else {
		trajectories <- list(trajectories=NULL)
		pop.observed <- get.pop.observed(pop.pred, country$code, sex=sex, age=age.idx)
		if(lage==max.age.idx) {
			if(sex == 'both') quant <- .get.pop.quantiles(pop.pred, '', adjust=adjust, ...)
			else quant <- .get.pop.quantiles(pop.pred, if(sex=='male') 'M' else 'F', adjust=adjust, ...)
		}
	}
	x1 <- names(pop.observed)[-length(pop.observed)]
	x2 <- if(!is.null(dimnames(trajectories$trajectories))) dimnames(trajectories$trajectories)[[1]]
			else dimnames(pop.pred$quantiles)[[3]]
	l <- length(x1) + length(x2)
	pred.table <- matrix(NA, ncol=2*length(pi)+1, nrow=l)
	rownames(pred.table) <- c(x1, x2)
	pred.table[x1,1] <- pop.observed[-length(pop.observed)]
	pred.table[x2,1] <- get.pop.traj.quantiles(quant, pop.pred, country$index, country$code, 
							trajectories=trajectories$trajectories,	q=0.5, sex=sex, age=age.idx, reload=is.null(expression), adjust=adjust)
	#if(is.na(pred.table[x2[1],1])) pred.table[x2[1],1] <- pop.observed[length(pop.observed)]
	
	colnames(pred.table) <- c('median', rep(NA,ncol(pred.table)-1))
	idx <- 2
	for (i in 1:length(pi)) {
		cqp <- get.pop.traj.quantiles(quant, pop.pred, country$index, country$code, 
							trajectories=trajectories$trajectories, pi=pi[i], sex=sex, age=age.idx, reload=is.null(expression), adjust=adjust)
		if (!is.null(cqp)) {
			pred.table[x2,idx:(idx+1)] <- t(cqp)
		} else{
			pred.table[x2,idx:(idx+1)] <- matrix(NA, nrow=l, ncol=2)
		}
		al <- (1-pi[i]/100)/2
		colnames(pred.table)[idx:(idx+1)] <- c(al, 1-al)
		idx <- idx+2
	}
	cn <- colnames(pred.table)[2:ncol(pred.table)]
	pred.table[,2:ncol(pred.table)] <- pred.table[,cn[order(cn)]]
	colnames(pred.table)[2:ncol(pred.table)] <- cn[order(cn)]
	if(half.child.variant && is.null(expression)) {
		# load the half child variants from trajectory file
		traj <- get.pop.trajectories(pop.pred, country$code, sex, age, nr.traj=0, adjust=adjust)
		if(!is.null(traj$half.child)) {
			pred.table <- cbind(pred.table, rbind(matrix(NA, nrow=length(x1), ncol=2), traj$half.child))
			colnames(pred.table)[(ncol(pred.table)-1):ncol(pred.table)] <- c('-0.5child', '+0.5child')
		}
	}
#	if(!is.null(adj.to.file)) {
#		pred.table[x2,] <- adjust.to.dataset(country$code, pred.table, adj.file=adj.to.file, years=x2, use='table')
#	}
	return(pred.table)
}

get.data.byage <- function(pop.pred, country.object, year=NULL, expression=NULL, pi=c(80, 95),
							sex=c('both', 'male', 'female'), nr.traj=NULL, typical.trajectory=FALSE, load.trajectories=TRUE) {
	projection <- TRUE
	projection.index <- 1
	if(!is.null(year)) {
		ind.proj <- get.predORobs.year.index(pop.pred, year)
		projection.index <- ind.proj['index']
		projection <- ind.proj['is.projection']
		if(is.null(projection.index)) stop('Year ', year, ' not found.')
	}
	country <- country.object
	cqp <- list()
	trajectories <- list(trajectories=NULL)
	last.open <- TRUE
	end.time.label <- FALSE
	if(projection) {
		if(!is.null(expression)) {
			trajectories <- get.pop.trajectories.from.expression.multiple.age(expression, pop.pred, nr.traj, 
										typical.trajectory=typical.trajectory)
			age.idx <- if(!is.null(rownames(trajectories$trajectories))) 
							as.integer(rownames(trajectories$trajectories)) 
						else 1:nrow(trajectories$trajectories)
		} else {
			if(load.trajectories) {
				trajectories <- get.pop.trajectories.multiple.age(pop.pred, country$code, sex, nr.traj=nr.traj, 
										typical.trajectory=typical.trajectory)
				age.idx <- trajectories$age.idx
			} else {
				if(sex == 'male') trajectories$quantiles <- pop.pred$quantilesMage
				else {
					if(sex=='female') trajectories$quantiles <- pop.pred$quantilesFage
					else # load trajectories to get quantiles 
						trajectories <- get.pop.trajectories.multiple.age(pop.pred, country$code, sex)
				}
				age.idx<-1:27
			}
			end.time.label <- TRUE
		}
		pop.median <- get.pop.traj.quantiles.byage(trajectories$quantiles, pop.pred, country$index, country$code, 
										projection.index, trajectories=trajectories$trajectories, q=0.5, sex=sex)
		ylim.loc <- c(min(pop.median, na.rm=TRUE), max(pop.median, na.rm=TRUE))
		if(!is.null(trajectories$trajectories)) 
			ylim.loc <- c(min(ylim.loc[1], trajectories$trajectories[,projection.index,trajectories$index], na.rm=TRUE),
							max(ylim.loc[2], trajectories$trajectories[,projection.index,trajectories$index], na.rm=TRUE))
		for (i in 1:length(pi))
			cqp[[i]] <- get.pop.traj.quantiles.byage(trajectories$quantiles, pop.pred, country$index, country$code, 
										projection.index, trajectories=trajectories$trajectories, pi=pi[i], sex=sex)
		if(length(pi) > 0)
			ylim.loc <- c(min(ylim.loc[1], sapply(cqp, min, na.rm=TRUE), na.rm=TRUE),
							max(ylim.loc[2], sapply(cqp, max, na.rm=TRUE), na.rm=TRUE))
		year.label <- get.pop.prediction.periods(pop.pred, end.time.only=end.time.label)[projection.index]
		if(length(age.idx)<27) last.open <- FALSE
	} else { # historical year 
		if(!is.null(expression)) {
			pop.observed <- get.pop.observed.from.expression.multiple.age(expression, pop.pred)
			pop.median <- pop.observed[,projection.index]
			age.idx <- if(!is.null(rownames(pop.observed))) as.integer(rownames(pop.observed)) else 1:nrow(pop.observed)
		} else {
			pop.observed <- get.pop.observed.with.age(pop.pred, country$code, sex=sex)
			pop.median <- pop.observed$data[,projection.index]
			age.idx <- pop.observed$age.idx
			end.time.label <- TRUE
		}		
		ylim.loc <- c(min(pop.median), max(pop.median))		
		year.label <- get.pop.observed.periods(pop.pred, end.time.only=end.time.label)[projection.index]
		if(length(age.idx)<21) last.open <- FALSE
	}
	age.labels <- get.age.labels(age.idx, age.is.index=TRUE, last.open=last.open)
	return(list(trajectories=trajectories, age.idx=age.idx, age.labels=age.labels, year.label=year.label, ylim.loc=ylim.loc,
					cqp=cqp, pop.median=pop.median, projection=projection, projection.index=projection.index))
}

pop.byage.plotAll <- function(pop.pred, 
									output.dir=file.path(getwd(), 'pop.byage'),
									output.type="png", expression=NULL, verbose=FALSE,  ...) {
	# plots pop trajectories by age for all countries
	if(!is.null(expression) && !grepl('XXX', expression, fixed=TRUE))
		stop('Expression must contain a mask "XXX" to be used as country code.')
	bayesTFR:::.do.plot.all.country.loop(pop.pred$countries[,'code'], meta=NULL, country.table=pop.pred$countries,
					output.dir=output.dir, func=.trajectories.plot.with.mask.replacement, output.type=output.type, 
					file.prefix='pop.byage.plot', plot.type='population graph', verbose=verbose, pop.pred=pop.pred, 
					expression=expression, plotfun='pop.byage.plot', ...)
}


pop.byage.plot <- function(pop.pred, country=NULL, year=NULL, expression=NULL, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), 
								  half.child.variant=FALSE,
								  nr.traj=NULL, typical.trajectory=FALSE,
								  xlim=NULL, ylim=NULL,  
								  xlab='', ylab='Population projection', main=NULL, 
								  lwd=c(2,2,2,1), col=c('red', 'red', 'blue', '#00000020'),
								  show.legend=TRUE, add=FALSE, ann=par('ann'), type='l', pch=NA, pt.cex=1, ...
								  ) {

	sex <- match.arg(sex)
	if (is.null(country)  && is.null(expression)) 
		stop('Argument "country" or "expression" must be given.')
		
	if(!is.null(country))
		country <- get.country.object(country, country.table=pop.pred$countries)
		
	data <- get.data.byage(pop.pred, country, year, expression, pi=pi,
							sex=sex, nr.traj=nr.traj, typical.trajectory=typical.trajectory)
	if(missing(ylab) && !is.null(expression)) ylab <- ''
	# recycle the last value to ensure there are at least 4 elements
	col <- c(col, rep(col[length(col)],3))
	lwd <- c(lwd, rep(lwd[length(lwd)],3))
	type <- c(type, rep(type[length(type)],3))
	pch <- c(pch, rep(pch[length(pch)],3))
	with(data, {
	x <- 1:length(age.idx)
	y <- pop.median
	if(is.null(xlim)) xlim <- c(min(x), max(x))
	if(is.null(ylim)) ylim <- ylim.loc
		
	if(is.null(main)) {
		if(!is.null(expression)) main <- paste(year.label, ':', expression)
		else {
			main <- paste(country$name, year.label) 
			if(sex != 'both') main <- paste(main, ': ', sex, sep='')		
		}
	}
	if(!add) {
		plot(x, y, type='n', xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
			panel.first = grid(), ann=ann,  xaxt='n', ...)
		if(ann) axis(1, at=1:length(age.idx), labels=age.labels, las=2)
	}
	# plot trajectories
	if(!is.null(trajectories$index) && !is.null(trajectories$trajectories)) {
		for (i in 1:length(trajectories$index)) {
			lines(x, trajectories$trajectories[,projection.index,trajectories$index[i]], 
					type=type[4], col=col[4], lwd=lwd[4], pch=pch[4], cex=pt.cex)
		}
	}
	# median
	lines(x,y, lwd=lwd[1], col=col[1], type=type[1], pch=pch[1], cex=pt.cex)
	# plot given CIs
	if(projection) {
		lty <- 2:(length(pi)+1)
		for (i in 1:length(pi)) {		
			if (!is.null(cqp[[i]])) {
				lines(x, cqp[[i]][1,], type=type[2], col=col[2], lty=lty[i], lwd=lwd[2], pch=pch[2], cex=pt.cex)
				lines(x, cqp[[i]][2,], type=type[2], col=col[2], lty=lty[i], lwd=lwd[2], pch=pch[2], cex=pt.cex)
			}
		}
		legend <- c('median', paste(pi, '% PI', sep=''))
	} else {
		legend <- c('observed')
		lty <- c()
		pi <- c()
	}
	lty <- c(1, lty)
	lwds <- c(lwd[1], rep(lwd[2], length(pi)))
	cols <- c(col[1], rep(col[2], length(pi)))
	types <- c(type[1], rep(type[2], length(pi)))
	pchs <- c(pch[1], rep(pch[2], length(pi)))
	if (half.child.variant && !is.null(trajectories$half.child)) {
		lty <- c(lty, max(lty)+1)
		llty <- length(lty)
		lines(x, trajectories$half.child[,projection.index,1], type=type[3], col=col[3], 
					lty=lty[llty], lwd=lwd[3], pch=pch[3], cex=pt.cex)
		lines(x, trajectories$half.child[,projection.index,2], type=type[3], col=col[3], 
					lty=lty[llty], lwd=lwd[3], pch=pch[3], cex=pt.cex)
		legend <- c(legend, '+/- 0.5 child')
		cols <- c(cols, col[3])
		lwds <- c(lwds, lwd[3])
		types <- c(types, type[3])
		pchs <- c(pchs, lwd[3])
	}
	if(show.legend && ann)
		legend('topleft', legend=legend, lty=lty, bty='n', col=cols, lwd=lwds, pch=pchs, pt.cex=pt.cex)
	})
}

pop.byage.table <- function(pop.pred, country=NULL, year=NULL, expression=NULL, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), 
								  half.child.variant=FALSE) {
	sex <- match.arg(sex)
	if (is.null(country)  && is.null(expression)) 
		stop('Argument "country" or "expression" must be given.')
		
	if(!is.null(country))
		country <- get.country.object(country, country.table=pop.pred$countries)
		
	data <- get.data.byage(pop.pred, country, year, expression, pi=pi,
							sex=sex, load.trajectories= if(half.child.variant && is.null(expression)) TRUE else FALSE)
	with(data, {
	x <- age.labels
	l <- length(x)
	pred.table <- matrix(NA, ncol=1, nrow=l)
	rownames(pred.table) <- x
	pred.table[x,1] <- pop.median
	colnames(pred.table) <- if(projection) paste('median', year.label, sep='.') else year.label
	if(projection && !is.null(cqp)) {
		idx <- 2
		for (i in 1:length(pi)) {
			pred.table <- cbind(pred.table, t(cqp[[i]]))
			al <- (1-pi[i]/100)/2
			colnames(pred.table)[idx:(idx+1)] <- c(al, 1-al)
			idx <- idx+2
		}
		cn <- colnames(pred.table)[2:ncol(pred.table)]
		pred.table[,2:ncol(pred.table)] <- pred.table[,cn[order(cn)]]
		colnames(pred.table)[2:ncol(pred.table)] <- cn[order(cn)]
	}
	if(projection && half.child.variant && is.null(expression) && !is.null(trajectories$half.child)) {
		# load the half child variants from trajectory file
		pred.table <- cbind(pred.table, trajectories$half.child[,projection.index,])
		colnames(pred.table)[(ncol(pred.table)-1):ncol(pred.table)] <- c('-0.5child', '+0.5child')
	}
	return(pred.table)
	})
}

"get.bPop.pyramid" <- function(data, ...) UseMethod("get.bPop.pyramid")

get.bPop.pyramid.matrix <- function(data, ...) get.bPop.pyramid.data.frame(data, ...)

get.bPop.pyramid.data.frame <- function(data, main.label=NULL, legend='observed', is.proportion=FALSE, ages=NULL, pop.max=NULL, 
										LRmain=c('Male', 'Female'), LRcolnames = c('male', 'female'), CI=NULL, ...) {
	if(is.null(colnames(data)) || any(!is.element(LRcolnames, colnames(data))))
		stop('Data must contain columns called ', paste(LRcolnames, collapse=', '), '.')
	if(!is.null(ages)) {
		if(length(ages) != nrow(data)) stop ("If ages is given, its length must be equal to the number of rows of data.")
		rownames(data) <- ages
	}
	if(is.na(is.proportion)) { # compute proportion
		data <- data/sum(data)
		is.proportion <- TRUE
	}
	maxx <- if(is.null(pop.max)) max(data) else pop.max
	pyr <- list(data)
	names(pyr) <- legend
	if(!is.null(CI)) {
		if(!is.list(CI))
			stop('CI must be a list with an element per confidence interval, each of which is a list of high and low data frames.')
		for(i in 1:length(CI)) {
			if(!all(is.element(c('high', 'low'), names(CI[[i]]))))
				stop('CI components must have elements high and low.')
			if(any(!is.element(LRcolnames, colnames(CI[[i]]$high))) || any(!is.element(LRcolnames, colnames(CI[[i]]$low))))
				stop('High and low CI must contain columns called ', paste(LRcolnames, collapse=', '), '.')
			if(is.null(pop.max)) maxx <- max(maxx, CI[[i]]$high, CI[[i]]$low)
		}

	}
	return(structure(list(
				label = main.label, 
				pyramid = pyr, CI = list(CI),
				is.proportion = is.proportion,
				pop.max=maxx,
				LRmain=LRmain,
				LRcolnames = LRcolnames), class='bayesPop.pyramid'))
}

get.bPop.pyramid.list <- function(data, main.label=NULL, legend=NULL, CI=NULL, ...) {
	if(is.null(legend)) legend <- names(data)
	ldata <- length(data)
	if(is.null(legend)) legend <- c('observed', rep('comparison', ldata-1))
	if(!is.null(CI) && length(CI) < ldata) {
		warning('CI must be the same length as data.')
		CI <- c(CI, vector('list', ldata-length(CI)))
	}
	pyr <- get.bPop.pyramid(data[[1]], main.label=main.label, legend=legend[[1]], ...)
	if(!is.null(data[[2]])) {
		for(i in 2:length(data)) {
			quant <- if(is.null(CI)) NULL else CI[[i]]
			pyri <- get.bPop.pyramid(data[[i]], legend=legend[[i]], CI=quant, ...)
			for(item in c('pyramid', 'CI'))
				pyr[[item]] <- c(pyr[[item]], pyri[[item]])
			pyr$pop.max <- max(pyr$pop.max, pyri$pop.max)
		}
	}
	return(pyr)
} 

.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}

get.bPop.pyramid.bayesPop.prediction <- function(data, country, year=NULL, indicator=c('P', 'B', 'D'),
												pi=c(80, 95), proportion=FALSE, age=1:21, 
												nr.traj=0, sort.pi=TRUE, pop.max=NULL, ...) {
	pop.pred <- data
	country <- get.country.object(country, country.table=pop.pred$countries)
	if(is.null(country$code)) stop('Country not found in the prediction object.')
	year.idx <- if(is.null(year)) 1 else sapply(year, get.prediction.year.index, pop.pred=pop.pred)
	lyears <- length(year.idx)
	draw.projection <- rep(TRUE, lyears)
	draw.projection[is.na(year.idx)] <- FALSE # if years not found in the prediction years, these are probably observed years
	draw.observed <- any(!draw.projection)
	pop.observed <- NULL
	indicator <- match.arg(indicator)
	if(draw.observed) {
		if(indicator == 'P')
			pop.observed <- get.pop.observed(pop.pred, country$code, sex='both')
		else
			pop.observed <- get.pop.observed.from.expression(paste(indicator, country$code, sep=''), pop.pred)
		year.idx[!draw.projection] <- sapply(year[!draw.projection], get.observed.year.index, pop.pred=pop.pred)
		if(all(is.na(year.idx))) stop('Unable to find data for year ', year)
		if(any(is.na(year.idx))) warning('Unable to find data for year ', year[is.na(year.idx)])
	}
	ages.idx <- age[age <=  length(pop.pred$ages)]
	lages <- length(ages.idx)
	nquant <- length(pi)
	quantiles.table <- trajs <- list(male=NULL, female=NULL)
	if(!any(draw.projection) || !any(draw.projection & (year.idx>1))) nquant <- 0
	if(nquant > 1 && sort.pi) pi<-sort(pi, decreasing=TRUE) # this is needed for drawing the largest intervals first (because of overlapping issues)
	ind.trajs <- NULL
	if(indicator != 'P') { # other indicators
		trajectoriesM <- get.popVE.trajectories.and.quantiles(pop.pred, country$code, 
								event=get.expression.indicators()[[indicator]], sex='male', 
								sum.over.ages=FALSE, q=c(0.5, (1-pi/100)/2, 1-(1-pi/100)/2))
		trajectoriesF <- get.popVE.trajectories.and.quantiles(pop.pred, country$code, 
								event=get.expression.indicators()[[indicator]], sex='female', 
								sum.over.ages=FALSE, q=c(0.5, (1-pi/100)/2, 1-(1-pi/100)/2))
		ages.idx <- ages.idx[is.element(ages.idx, trajectoriesM$age.idx.raw)]
		ages.idx.q <- trajectoriesM$age.idx[is.element(trajectoriesM$age.idx.raw, ages.idx)]
		lages <- length(ages.idx)
		ind.trajs <- list(male=trajectoriesM, female=trajectoriesF)
		for(sex in c("male", "female")) # drop ages that are not needed
			ind.trajs[[sex]]$trajectories <- ind.trajs[[sex]]$trajectories[ages.idx.q,,,drop=FALSE]
		if(proportion) {
			sumTrajs <- apply(trajectoriesM$trajectories + trajectoriesF$trajectories, c(2,3), sum)
			if(is.null(dim(sumTrajs))) sumTrajs <- abind(sumTrajs, along=2)						
			for (itraj in 1:dim(trajectoriesM$trajectories)[3]) {
				sTm <- matrix(rep(sumTrajs[,itraj], lages), nrow=lages, byrow=TRUE)
				ind.trajs[['male']]$trajectories[,,itraj] <- ind.trajs[['male']]$trajectories[,,itraj]/sTm
				ind.trajs[['female']]$trajectories[,,itraj] <- ind.trajs[['female']]$trajectories[,,itraj]/sTm
			}
		} else 
			quantiles.table <- list(male=trajectoriesM$quantiles, female=trajectoriesF$quantiles)

		
	} else # population indicator
		quantiles.table <- if(proportion) list(male=pop.pred$quantilesPropMage, female=pop.pred$quantilesPropFage)
                       else list(male=pop.pred$quantilesMage, female=pop.pred$quantilesFage)
    is.valid.pi <- if(proportion && nquant>0) is.saved.pi(pop.pred, pi)
                   else rep(TRUE, nquant)
	maxx<-0
	age.labels <- get.age.labels(pop.pred$ages[ages.idx])
	true.nquant <- sum(is.valid.pi)
	pyr <- list()
	pyr.ci <- list()
	for(yi in 1:lyears) {
		pyr[[yi]] <- data.frame(male=rep(NA, lages), female=rep(NA, lages), row.names=age.labels)
		pyr.ci[[yi]] <- list()
		if(!draw.projection[yi] || (draw.projection[yi] && year.idx[yi] == 1)) next
		pyr.ci[[yi]] <- vector("list", length=true.nquant)
		if(true.nquant > 0) {
			names(pyr.ci[[yi]]) <- pi[is.valid.pi]
			for(tpi in names(pyr.ci[[yi]])) 
				pyr.ci[[yi]][[tpi]] <- list(low=data.frame(male=rep(NA, lages), female=rep(NA, lages), row.names=age.labels),
										high=data.frame(male=rep(NA, lages), female=rep(NA, lages), row.names=age.labels))
		}
	}
	for(sex in c('male', 'female')) {
		dimt <- dim(quantiles.table[[sex]])
		dimn <- dimnames(quantiles.table[[sex]])
		this.trajs <- table <- NULL
		for(iage in 1:lages) {
			if(any(draw.projection)) {
				if(length(dim(quantiles.table[[sex]]))==4) { # population
					table <- drop(quantiles.table[[sex]][,ages.idx[iage],,])
					table <- array(table, dimt[c(1,3:4)], dimnames=c(list(NULL), dimn[3], dimn[4]))
					cidx <- country$index
					ci.reload <- TRUE
				} else { # other indicator
					if(!is.null(quantiles.table[[sex]])) {
						table <- drop(quantiles.table[[sex]][ages.idx.q[iage],,])
						table <- array(table, dimt[c(2,3)], dimnames=c(dimn[2], dimn[3]))
					} else this.trajs <- ind.trajs[[sex]]$trajectories[ages.idx.q[iage],,]
					cidx <- NULL
					ci.reload <- FALSE
				}
				med <- get.pop.traj.quantiles(table, pop.pred, cidx, country$code, trajectories=this.trajs,
												q=0.5, reload=FALSE, sex=sex, age=ages.idx[iage])
			}
			if(any(!draw.projection)) {
				if(indicator != 'P') {
					observed.data <- get.pop.observed.from.expression(
							paste(indicator, country$code, '_', substr(toupper(sex), 1,1), '[', ages.idx[iage], ']', sep=''), pop.pred)
				} else observed.data <- get.pop.observed(pop.pred, country$code, sex=sex, age=ages.idx[iage])
			}
			for(yi in 1:lyears) {				
				pyr[[yi]][iage,sex] <- if(draw.projection[yi]) med[year.idx[yi]] 
											else observed.data[year.idx[yi]]/(if(proportion) pop.observed[year.idx[yi]] else 1)
				if(is.na(pyr[[yi]][iage,sex])) pyr[[yi]][iage,sex] <- 0
				maxx <- max(maxx, pyr[[yi]][iage,sex])
			}
			if (nquant == 0) next
			for (i in 1:nquant) {
				if (!is.valid.pi[i]) next
				pi.name <- as.character(pi[i])
				quant <- get.pop.traj.quantiles(table, 
												pop.pred, cidx, country$code, trajectories=this.trajs,
												pi=pi[i], reload=ci.reload, sex=sex, age=ages.idx[iage])
				for(yi in 1:lyears) {
					if(draw.projection[yi] & (year.idx[yi] > 1)) {
						pyr.ci[[yi]][[pi.name]]$low[iage,sex] <- quant[1,year.idx[yi]]
						pyr.ci[[yi]][[pi.name]]$high[iage,sex] <- quant[2,year.idx[yi]]
						maxx <- max(maxx, pyr.ci[[yi]][[pi.name]]$high[iage,sex], pyr.ci[[yi]][[pi.name]]$low[iage,sex])
					}
				}
			}
		}
	}
	trajs <- list()
	male.trajectories <- female.trajectories <- NULL
	if((is.null(nr.traj) || nr.traj > 0) && any(draw.projection & (year.idx > 1))) {
		mtraj <- ind.trajs$male
		if(is.null(mtraj)) 
			mtraj <- get.pop.trajectories.multiple.age(pop.pred, country$code, sex='male', 
										age=ages.idx, nr.traj, proportion=proportion)
		if(!is.null(dim(mtraj$trajectories))) {
			ftraj <- ind.trajs$female
			if(is.null(ftraj)) 
				ftraj <- get.pop.trajectories.multiple.age(pop.pred, country$code, sex='female', 
										age=ages.idx, nr.traj, proportion=proportion) 
			for(yi in 1:lyears) {
				if(!draw.projection[yi] || (year.idx[yi] == 1)) {
					trajs[[yi]] <- list()
					next
				}
				male.trajectories <- drop(mtraj$trajectories[,year.idx[yi],mtraj$index])
				male.trajectories <- array(male.trajectories, c(dim(mtraj$trajectories)[1],length(mtraj$index)), 
												dimnames=list(age.labels, NULL))
				female.trajectories <- drop(ftraj$trajectories[,year.idx[yi],ftraj$index])
				female.trajectories <- array(female.trajectories, c(dim(ftraj$trajectories)[1],length(ftraj$index)),
												dimnames=list(age.labels, NULL))
				if(!is.null(male.trajectories)) maxx <- max(maxx, male.trajectories, female.trajectories)
				trajs[[yi]] <- list(male=male.trajectories, female=female.trajectories)
			}
		}
	}
	names(pyr)[1] <- if(draw.projection[1] && year.idx[1] > 1) 'median' else 'observed'
	if(indicator=='P') {
		proj.years <- litem('proj.years.pop', pop.pred, pop.pred$proj.years+2)
		obs.years <- as.integer(names(pop.observed))
		main.year <- (if(draw.projection[1]) proj.years else obs.years)[year.idx[1]]
		if(lyears > 1)
			names(pyr)[2:lyears] <- ifelse(draw.projection[2:lyears], proj.years[year.idx[2:lyears]],  
											obs.years[year.idx[2:lyears]])
	} else {
		proj.years <- pop.pred$proj.years
		obs.years <- as.integer(names(pop.observed))
		main.year <- paste((if(draw.projection[1]) proj.years else obs.years)[year.idx[1]] + c(-3, 2), collapse='-')
		if(lyears > 1)
			names(pyr)[2:lyears] <- sapply(lapply(ifelse(draw.projection[2:lyears], proj.years[year.idx[2:lyears]],  
											obs.years[year.idx[2:lyears]]), '+', c(-3, 2)), paste, collapse='-')
	}
	indicator.name <- if(indicator=='P') 'Population' else .simpleCap(get.expression.indicators()[[indicator]])
			
	return(structure(list(
				label = paste(indicator.name, ' in ', country$name, ': ', main.year, sep=''), 
				pyramid = pyr, CI = pyr.ci,
				trajectories = if(length(trajs) > 0) trajs else NULL,
				is.proportion = proportion,
				pop.max=if(is.null(pop.max)) maxx else pop.max,
				LRmain=c('Male', 'Female'),
				LRcolnames = c('male', 'female')
				), class='bayesPop.pyramid'))
}


plot.bayesPop.pyramid <- function(x, ...) {
	if(is.null(x$trajectories))
		pop.pyramid(x, ...)
	else pop.trajectories.pyramid(x, ...)
}

"pop.pyramid" <- function(pop.object, ...) UseMethod("pop.pyramid")

pop.pyramid.bayesPop.pyramid <- function(pop.object, main=NULL, show.legend=TRUE, 
										pyr1.par=list(border='black', col=NA, density=NULL, height=0.9),
										pyr2.par =list(density=-1, height=0.3), 
										col.pi = NULL, ann=par('ann'), axes=TRUE, grid=TRUE, 
										cex.main=0.9, cex.sub=1, cex=1, cex.axis=1, ...) {
	mgp <- par('mgp')
	mar <- par('mar')
	par(mgp=c(3,0.5,0)) 
	par(mar=c(5, 4, 2, 4) + 0.1)
	if((is.null(pop.object$pyramid) || length(pop.object$pyramid) == 0) && is.null(pop.object$CI)) 
		stop('Nothing to be plotted. Either pyramid or CI must be given in pop.object.')
	age.labels <- rownames(if(!is.null(pop.object$pyramid[[1]])) pop.object$pyramid[[1]] 
						   else {if(!is.null(pop.object$pyramid[[2]])) pop.object$pyramid[[2]] else pop.object$CI[[1]][[1]]$low})
	if(is.null(age.labels)) stop('Row names must be given to determine age labels.')
	lages <- length(age.labels)
	nquant <- length(pop.object$CI[[1]])
	draw.past <- (length(pop.object$pyramid) > 1) && !is.null(pop.object$pyramid[[2]])
	draw.median <- !is.null(pop.object$pyramid[[1]])
	pyr1.par.default <- list(border='black', col=NA, density=NULL, height=0.9)
	for(item in names(pyr1.par.default)) if(is.null(pyr1.par[[item]])) pyr1.par[[item]] <- pyr1.par.default[[item]]
	col.pyr2.default <- adjustcolor('#f46d43', alpha.f=0.3)
	pyr2.par.default <- list(border=col.pyr2.default, col=col.pyr2.default, density=-1, height=0.3)
	for(item in names(pyr2.par.default)) if(is.null(pyr2.par[[item]])) pyr2.par[[item]] <- pyr2.par.default[[item]]
	pyr1.half.height <- pyr1.par[['height']]/2
	pyr2.half.height <- pyr2.par[['height']]/2
	pyr1q.half.height <- max(0.1, min(1, pyr1.half.height + pyr1.half.height/4))
	pyr2q.half.height <- max(0.1, min(1, pyr2.half.height + pyr2.half.height/4))
	cols <- lwd <- legend <- c()
	quantiles <- pop.object$CI[[1]]
	pyr1 <- pop.object$pyramid[[1]]
	pyr2 <- if(draw.past) pop.object$pyramid[[2]] else NULL
	#if(draw.median)
	#	cohort.labels <- if("_cohorts_" %in% colnames(pyr1)) pyr1[,"_cohorts_"] else age.labels
	#else cohort.labels <- if(!is.null(pyr2) && "_cohorts_" %in% colnames(pyr2)) pyr2[,"_cohorts_"] else age.labels
	with(pop.object, {
		maxx <- pop.max	
		proportion <- !is.null(is.proportion) && is.proportion
		male <- LRcolnames[1]
		female <- LRcolnames[2]
		plot(c(-maxx, maxx), c(-0.5, lages-0.5), type='n', axes=FALSE, xlab = "", ylab = "", ann=ann)
		if(ann) mtext(LRmain, at=c(-maxx/2, maxx/2), side=3, cex=cex.sub)
		age.axis.at <- 0:(lages-1)
		labels <- .get.xtick.labels.for.pyramid(maxx, proportion) 
		xat <- c(-labels, labels[1:length(labels)])
		if(axes) {
			axis(1, at=xat, labels=c(labels, labels[1:length(labels)]), cex.axis=cex.axis)
			axis(2, at=age.axis.at, labels=age.labels, las=2, cex.axis=cex.axis)
			axis(4, at=age.axis.at, labels=age.labels, las=2, cex.axis=cex.axis)
		}
		if(grid) {#grid(length(labels))
			gridxat <- xat[seq(1, length(xat), by=2)]
			segments(gridxat, rep(-0.5-lages/25, length(gridxat)), gridxat, lages-1+lages/25, col="lightgray", lty = "dotted")
			gridyat <- age.axis.at[seq(1, lages, by=2)]
			segments(rep(-maxx-2*maxx/25, length(gridyat)), gridyat, rep(maxx+2*maxx/25, length(gridyat)), 
						gridyat, col="lightgray", lty = "dotted")
		}
		if(nquant > 0) {			
			if(is.null(col.pi)) {
				if(nquant < 10) {
					# from RcolorBrewer: brewer.pal(9, "YlGnBu")
					col.pi.default <- c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
					if(nquant < 5) # remove the extreme colors at both ends
						col.pi.default <- col.pi.default[2:(length(col.pi.default)-3)]
					cols <- if(nquant == 2) c("#C7E9B4", "#41B6C4") else col.pi.default[sort(seq(length(col.pi.default), 1, length=nquant))]
				} else { # more than 9 PIs
					cols <- rainbow(nquant, start=0.15)
				}
			} else col <- rep(col.pi, nquant)[1:nquant]
			for(i in 1:nquant) {
				rect(-quantiles[[i]]$high[,male], age.axis.at-pyr1q.half.height, 
						-quantiles[[i]]$low[,male], age.axis.at+pyr1q.half.height, col=cols[i],
						border= NA)
				rect(quantiles[[i]]$low[,female], age.axis.at-pyr1q.half.height, 
						quantiles[[i]]$high[,female], 
				 		age.axis.at+pyr1q.half.height, col=cols[i], border=NA)
			}
			legend <- c(legend, paste(names(quantiles), '% PI', sep=''))
			lwd <- c(lwd, rep(5, nquant))
		}
		if(draw.median) {
			rect(-pyr1[,male], age.axis.at-pyr1.half.height, rep(0, lages), age.axis.at+pyr1.half.height,
					col=pyr1.par$col, border=pyr1.par$border, density=pyr1.par$density)
			segments(-pyr1[,male], age.axis.at-pyr1.half.height, -pyr1[,male], age.axis.at+pyr1.half.height, 
					col=pyr1.par$border, lwd=3)
			rect(rep(0, lages), age.axis.at-pyr1.half.height, pyr1[,female], age.axis.at+pyr1.half.height, #lwd=2,
					col=pyr1.par$col, border=pyr1.par$border, density=pyr1.par$density)
			segments(pyr1[,female], age.axis.at-pyr1.half.height, pyr1[,female], age.axis.at+pyr1.half.height, 
					col=pyr1.par$border, lwd=3)
			lwd <- c(3, lwd)
			cols <- c(pyr1.par$border, cols)
			legend <- c(names(pyramid)[1], legend)
		}		
		if(draw.past) {
			rect(-pyr2[,male], age.axis.at-pyr2.half.height, rep(0, lages), age.axis.at+pyr2.half.height,
					col=pyr2.par$col, border=pyr2.par$border, density=pyr2.par$density)
			segments(-pyr2[,male], age.axis.at-pyr2.half.height, -pyr2[,male], age.axis.at+pyr2.half.height, 
					col=pyr2.par$border, lwd=3)
			rect(rep(0, lages), age.axis.at-pyr2.half.height, pyr2[,female], age.axis.at+pyr2.half.height,
					col=pyr2.par$col, border=pyr2.par$border, density=pyr2.par$density)
			segments(pyr2[,female], age.axis.at-pyr2.half.height, pyr2[,female], age.axis.at+pyr2.half.height, 
					col=pyr2.par$border, lwd=3)
			legend <- c(legend, names(pyramid)[2])
    		cols <- c(cols, pyr2.par$border)
    		lwd <- c(lwd, 3)
		}
		lines(c(0,0), c(age.axis.at[1]-pyr1.half.height, age.axis.at[length(age.axis.at)]+pyr1.half.height), col='black')	
		if(show.legend && ann) legend('topright', legend=legend, bty='n', col=cols, lwd=lwd, cex=cex)
		if(is.null(main)) main <- if(exists('label')) label else ""
		if(ann) title(main, line=1, cex.main=cex.main)
	})	
	par(mgp=mgp, mar=mar)
}

pop.pyramid.bayesPop.prediction <- function(pop.object, country, year=NULL, indicator=c('P', 'B', 'D'),
											pi=c(80, 95), proportion=FALSE,
											age=1:21, plot=TRUE, pop.max=NULL, ...) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	data <- get.bPop.pyramid(pop.object, country, year=year, indicator=indicator, pi=pi, proportion=proportion, age=age, pop.max=pop.max)
	if (plot) pop.pyramid(data, ...)
	invisible(data)
}

.do.plot.pyramid.all <- function(pop.pred, output.dir, func, year=NULL, output.type="png", 
						file.prefix='pyr', plot.type='pyramid(s)', one.file=FALSE,
						main=NULL, verbose=FALSE, ...) {
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.countries <- pop.pred$countries[,'name']
	if(one.file && !(output.type %in% c("pdf", 'postscript'))) output.type <- "pdf"
	postfix <- output.type	
	if(output.type=='postscript') postfix <- 'ps'
	if(is.null(year)) year <- list(pop.pred$present.year)
	if(!is.list(year)) year <- list(year)
	main.arg <- main
	if(one.file) 
		do.call(output.type, list(file.path(output.dir, 
										paste(file.prefix, '.', postfix, sep=''))))
	for (country in all.countries) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries)
		if(verbose)
			cat('Creating', plot.type, 'for', country, '(', country.obj$code, ')\n')
		if(!is.null(main) && grepl('XXX', main, fixed=TRUE))
			main.arg <- gsub('XXX', as.character(country.obj$name), main, fixed=TRUE)
		for(y in year) {
			if(!one.file) do.call(output.type, list(file.path(output.dir, 
										paste(file.prefix, paste(y, collapse='_'), '_c', country.obj$code, '.', postfix, sep=''))))
			do.call(func, list(pop.pred, country=country.obj$code, year=y, main=main.arg, ...))
			if(!one.file) dev.off()
		}
	}
	if(one.file) dev.off()
	if(verbose)
		cat('\nPyramids stored into', output.dir, '\n')
}


pop.pyramidAll <- function(pop.pred, year=NULL,
									output.dir=file.path(getwd(), 'pop.pyramid'),
									output.type="png", one.file=FALSE, verbose=FALSE, ...) {
	# plots pyramid for all countries
	.do.plot.pyramid.all(pop.pred, output.dir, pop.pyramid, year=year, output.type=output.type, 
				one.file=one.file, verbose=verbose, ...)
}


"pop.trajectories.pyramid" <- function(pop.object, ...) UseMethod("pop.trajectories.pyramid")

pop.trajectories.pyramid.bayesPop.prediction <- function(pop.object, country, year=NULL, indicator=c('P', 'B', 'D'), 
														pi=c(80, 95), nr.traj=NULL, proportion=FALSE, 
														age=1:21, plot=TRUE, pop.max=NULL, ...) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	data <- get.bPop.pyramid(pop.object, country, year=year, indicator=indicator, pi=pi, nr.traj=nr.traj, proportion=proportion, 
							age=age, sort.pi=FALSE, pop.max=pop.max)
	if(plot) pop.trajectories.pyramid(data, ...)
	invisible(data)
}

pop.trajectories.pyramid.bayesPop.pyramid  <- function(pop.object, main=NULL, show.legend=TRUE, 
													col=rainbow, col.traj='#00000020',
													lwd=2, ann=par('ann'), axes=TRUE, grid=TRUE, 
													cex.main=0.9, cex.sub=1, cex=1, cex.axis=1, ...) {
	# col/lwd is color and line width for:
	# 1. median, 2. quantiles, 3. past data, 4. trajectories
	mgp <- par('mgp')
	oma <- par('oma')
	mar <- par('mar')
	#par(mfrow=c(1,2),  mar=c(5,6,2,-0.1)+0.1)
	par(oma = c(0, 0, 2, 0), mgp=c(3,0.5,0))
	par(mar=c(5, 4, 2, 4) + 0.1)
	if((is.null(pop.object$pyramid) || length(pop.object$pyramid) == 0) && is.null(pop.object$CI) && is.null(pop.object$trajectories))
		stop('Nothing to be plotted. Either pyramid, CI or trajectories must be given in pop.object.')
	pyr.indicator <- !sapply(pop.object$pyramid, is.null)
	lpyr <- length(pyr.indicator)
	ci.indicator <- if(is.null(pop.object$CI)) rep(FALSE, lpyr) 
						else sapply(pop.object$CI, function(x) return(length(x)>0))
	traj.indicator <- if(is.null(pop.object$trajectories)) rep(FALSE, lpyr) 
						else sapply(pop.object$trajectories, function(x) return(length(x)>0))
	age.labels <- rownames(if(any(pyr.indicator)) pop.object$pyramid[[which(pyr.indicator)[1]]]
						   else {if(any(ci.indicator)) pop.object$CI[[which(ci.indicator)[1]]]$low else {
						   			if(any(traj.indicator)) pop.object$trajectories[[which(traj.indicator)[1]]]$male else NULL}})
	if(is.null(age.labels)) stop('Row names must be given to determine age labels.')
	lages <- length(age.labels)
	if(length(lwd) < lpyr) {
		lwd <- rep(lwd, lpyr)
		lwd <- lwd[1:lpyr]
	}
	if(is.function(col)) col <- do.call(col, list(lpyr))
	else {
		if(length(col) < lpyr) {
			col <- rep(col, lpyr)
			col <- col[1:lpyr]
		}
	}
	if(length(col.traj) < lpyr) {
		col.traj <- rep(col.traj, lpyr)
		col.traj <- col.traj[1:lpyr]
	}
	legend <- lty <- cols <- lwds <- ltys <- c()
	maxx <- pop.object$pop.max
	proportion <- !is.null(pop.object$is.proportion) && pop.object$is.proportion
	labels <- .get.xtick.labels.for.pyramid(maxx, proportion)
	age.axis.at <- 0:(lages-1)
	plot(c(-maxx,maxx), c(0, lages-0.5), type='n', axes=FALSE, xlab = "", ylab = "", ann=ann)
	if(ann) mtext(pop.object$LRmain, at=c(-maxx/2, maxx/2), side=3, cex=cex.sub)
	xat <- c(-labels, labels[1:length(labels)])
	if(axes) {
		axis(1, at=xat, labels=c(labels, labels[1:length(labels)]), cex.axis=cex.axis)
		axis(2, at=age.axis.at, labels=age.labels, las=2, cex.axis=cex.axis)
		axis(4, at=age.axis.at, labels=age.labels, las=2, cex.axis=cex.axis)
	}
	if(grid) {#grid(length(labels))
		gridxat <- xat[seq(1, length(xat), by=2)]
		segments(gridxat, rep(-0.5-lages/25, length(gridxat)), gridxat, lages-1+lages/25, col="lightgray", lty = "dotted")
		gridyat <- age.axis.at[seq(1, lages, by=2)]
		segments(rep(-maxx-2*maxx/25, length(gridyat)), gridyat, rep(maxx+2*maxx/25, length(gridyat)), 
						gridyat, col="lightgray", lty = "dotted")
	}
	lines(c(0,0), c(-lages/25, lages))

	with(pop.object, {
		for(ipyr in 1:length(pyr.indicator)) {
			nquant <- if (ci.indicator[ipyr]) length(CI[[ipyr]]) else 0
			draw.median <- pyr.indicator[ipyr]
			draw.traj <- traj.indicator[ipyr]
			male <- LRcolnames[1]
			female <- LRcolnames[2]
			nr.traj <- if(draw.traj) ncol(pop.object$trajectories[[ipyr]][[male]]) else 0				
			if(draw.traj) 
				for(i in 1:nr.traj) {
					lines(-trajectories[[ipyr]][[male]][,i], age.axis.at, col=col.traj[ipyr], lwd=1)
					lines(trajectories[[ipyr]][[female]][,i], age.axis.at, col=col.traj[ipyr], lwd=1)
				}
			if(draw.median) {
				lines(-pyramid[[ipyr]][,male], age.axis.at, col=col[ipyr], lwd=lwd[ipyr])
				lines(pyramid[[ipyr]][,female], age.axis.at, col=col[ipyr], lwd=lwd[ipyr])
				cols <- c(cols, col[ipyr])
				lwds <- c(lwds, lwd[ipyr])
				ltys <- c(ltys, 1)
				legend <- c(legend, names(pyramid)[ipyr])
			}
			if(nquant > 0) {
				lty <- 2:(nquant+1)
				for(i in 1:nquant) {
					lines(-CI[[ipyr]][[i]]$low[,male], age.axis.at, col=col[ipyr], lwd=lwd[ipyr], lty=lty[i])
					lines(-CI[[ipyr]][[i]]$high[,male], age.axis.at, col=col[ipyr], lwd=lwd[ipyr], lty=lty[i])
					lines(CI[[ipyr]][[i]]$low[,female], age.axis.at, col=col[ipyr], lwd=lwd[ipyr], lty=lty[i])
					lines(CI[[ipyr]][[i]]$high[,female], age.axis.at, col=col[ipyr], lwd=lwd[ipyr], lty=lty[i])
				}
				cols <- c(cols, rep(col[ipyr], nquant))
				lwds <- c(lwds, rep(lwd[ipyr], nquant))
				legend <- c(legend, paste(names(CI[[ipyr]]), '% PI', sep=''))
				ltys <- c(ltys, lty)
			}
		}
		if(is.null(main)) main <- if(exists("label")) label else ""
		if(ann) title(main, cex.main=cex.main, line=1)	
		if(show.legend && ann) legend('topright', legend=legend, lty=ltys, bty='n', col=cols, lwd=lwds, cex=cex)
	})
	par(mgp=mgp, oma=oma, mar=mar)
}

pop.trajectories.pyramidAll <- function(pop.pred, year=NULL,
									output.dir=file.path(getwd(), 'pop.traj.pyramid'),
									output.type="png", one.file=FALSE, verbose=FALSE, ...) {
	# plots pyramid for all countries and all years given by 'year'
	.do.plot.pyramid.all(pop.pred, output.dir, pop.trajectories.pyramid, year=year, output.type=output.type, 
					one.file=one.file, plot.type='trajectory pyramid(s)', verbose=verbose, ...)
}


.get.xtick.labels.for.pyramid <- function(maxx, proportion){
	if(proportion) {
		rmaxx <- round(maxx*100)
		nticks <- min(rmaxx+1,11)
		dec <- 2
		if(rmaxx <= 5) {
			nticks <- 2*rmaxx + 1
			dec <- 3
		}
		return(round(seq(0, maxx, length=nticks), dec))
	}
	return(round(signif(seq(0, maxx, length=min(11,maxx)),2)))
}


get.data.for.worldmap.bayesPop.prediction <- function(pred, quantile=0.5, year=NULL,
							projection.index=1, pi=NULL, sex=c('both', 'male', 'female'), age='all', expression=NULL, ...) {
	quantiles <- quantile
	if (!is.null(pi)) {
		qlower <- (1-pi/100)/2
		quantiles <- c(quantile, qlower, 1-qlower)
	}
	if(!all(is.element(as.character(quantiles), as.character(get.quantiles.to.keep()))))
		stop('The quantile and pi arguments must correspond to the following quantiles:\n', get.quantiles.to.keep())
	sex <- match.arg(sex)
	projection <- TRUE
	if(!is.null(year)) {
		ind.proj <- get.predORobs.year.index(pred, year)
		projection.index <- ind.proj['index']
		projection <- ind.proj['is.projection']
		if(is.null(projection.index)) stop('Projection year ', year, ' not found.')
	}
	if(projection) {
		if(!is.null(expression)) {
			data <- get.pop.from.expression.all.countries(expression, pred, quantiles, projection.index)
		} else {
			if(!all(is.element(as.character(quantiles), dimnames(pred$quantiles)[[2]])))
				stop('Some of the quantiles ', paste(quantiles, collapse=', '), ' not found.\nAvailable: ', 
							paste(dimnames(pred$quantiles)[[2]], collapse=', '), 
					 '\nCheck arguments "quantile" and "pi".')
			data <- get.pop.all.countries(pred, quantiles, projection.index, sex, age)
			if(is.null(data)) stop('Maps are not supported for the given combination of sex and age.')
		}
		period <- get.pop.prediction.periods(pred)[projection.index]
	} else { # observed data
		data <- if(!is.null(expression))
					get.pop.observed.from.expression.all.countries(expression, pred, projection.index)
				else 
					get.pop.observed.all.countries(pred, projection.index, sex=sex, age=age)
		period <- get.pop.observed.periods(pred)[projection.index]
	}
	rownames(data) <- NULL
	low<-NULL
	up<-NULL
	res <- data
	if(!is.null(dim(data))) {		
		res <- data[,1]
		if(ncol(data) > 1) {
			low <- data[,2]
			up <- data[,3]
		}
	}
	return(list(period=period, data=res, country.codes=pred$countries$code, lower=low, upper=up, 
				sex=sex, age=age, expression=expression))
}


pop.map <- function(pred, sex=c('both', 'male', 'female'), age='all', expression=NULL, ...) {
	if(!requireNamespace("rworldmap", quietly=TRUE)) {
		warning("Package 'rworldmap' is not installed. If 'googleVis' is installed, use pop.map.gvis(...).")
		return()
	}
	return(bayesTFR::tfr.map(pred, par.name=expression, data.args=list(sex=sex, age=age, expression=expression), ...))
}

.map.main.default.bayesPop.prediction <- function(pred, dp, ...) {
	if(!is.null(dp$expression)) main <- dp$expression
	else {
		main <- 'Pop'
		if(dp$sex != 'both') main <- paste(main, dp$sex)
		if(dp$age[1] != 'all') main <- paste(main, get.age.labels(dp$age, collapsed=TRUE, age.is.index=TRUE))
	}
	return(paste(main,': quantile', sep=''))
}

get.pop.map.parameters <- function(pred, expression=NULL, sex=c('both', 'male', 'female'), age='all', 
								range=NULL, nr.cats=50, same.scale=TRUE, quantile=0.5, ...) {
	map.pars <- list(pred=pred, quantile=quantile, expression=expression, sex=sex, age=age, ...)
	if (same.scale) {
		if(is.null(range)) {
			if(!all(is.element(quantile, get.quantiles.to.keep()))) 
			stop('The quantile and pi arguments must correspond to the following quantiles:\n', get.quantiles.to.keep())
			sex <- match.arg(sex)
			data <- if(!is.null(expression)) get.pop.from.expression.all.countries(expression, pred, quantile, 1)
					else get.pop.all.countries(pred, quantile, 1, sex=sex, age=age)	
			range <- c(min(data), max(data))
		}
		quant.values <- seq(range[1], range[2], length=nr.cats)
		col <- rainbow(500, start=0, end=0.67)[seq(500, 1, length=length(quant.values)-1)]
		map.pars$catMethod <- quant.values
	} else {
		col <- rainbow(500, start=0, end=0.67)[seq(500, 1, length=nr.cats)]
		map.pars$numCats <- nr.cats
	}
	map.pars$colourPalette <- col
	return(map.pars)		
}

pop.map.gvis <- function(pred, ...){
	if(!requireNamespace("googleVis", quietly=TRUE)) {
		warning("Package 'googleVis' is not installed. If 'rworldmap' is installed, use pop.map(...).")
		return()
	}
	bdem.map.gvis(pred, ...)
}
					
bdem.map.gvis.bayesPop.prediction <- function(pred,  ...) {
	bayesTFR:::.do.gvis.bdem.map('pop', 'Population', pred, ...)
}

pop.cohorts.plot <- function(pop.pred, country=NULL, expression=NULL, cohorts=NULL, cohort.data=NULL, pi=c(80, 95), 
								dev.ncol=5, show.legend=TRUE, legend.pos="bottomleft", ann=par('ann'), add=FALSE, 
								xlab="", ylab="", main=NULL, xlim=NULL, ylim=NULL, col='red', ...) {
	.round.to.lower5 <- function(x) 5*floor(x/5) 

	if(is.null(cohort.data)) 
		cohort.data <- cohorts(pop.pred, country=country, expression=expression, pi=pi)
	all.cohorts <- names(cohort.data)[-which(names(cohort.data) == 'last.observed')]
	all.cohorts.num.start <- as.integer(substr(all.cohorts, 1, 4))
	if(is.null(cohorts)) 
		cohorts <- seq(cohort.data[['last.observed']], by=5, 
							length=min(10, sum(all.cohorts.num.start > cohort.data[['last.observed']])))
	if(any(is.numeric(cohorts))) {
		# convert to the from-to format, e.g. 2000-2005
		from.cohorts <- .round.to.lower5(cohorts)
		cohorts <- paste(from.cohorts, '-', from.cohorts+5, sep="")
	}
	if(is.null(xlim))
		xlim <- range(unlist(sapply(cohorts, function(x) range(as.integer(colnames(cohort.data[[x]]))))))
	if(is.null(ylim))
		ylim <- range(unlist(sapply(cohorts, function(x) range(cohort.data[[x]]))))

	nplots <- length(cohorts)
	if (nplots < dev.ncol) {
       	ncols <- nplots
		nrows <- 1
	} else {
		ncols <- dev.ncol
		nrows <- ceiling(nplots/dev.ncol)
	}

	pi <- sort(pi)
	legend <- c('median', paste(pi, '% PI', sep=''))
	lty <- c(1, 2:(length(pi)+1))
	cols <- rep(col, 1+length(pi))
	if(nplots > 1) {
		cur.mgp <- par('mgp')
		cur.oma <- par('oma')
		cur.mar <- par('mar')
		par(mfrow=c(nrows,ncols),  oma = c(0, 0, 2, 0))
		par(mar=c(2,2,1,0.4)+0.1, mgp=c(1,0.3,0))
	}
	for(iplot in 1:nplots) {
		this.data <- cohort.data[[cohorts[iplot]]]
		this.x <- as.integer(colnames(this.data))
		if(!add)
			plot(this.x, this.data["median",], type="l", col=col, xlim=xlim, ylim=ylim,
				ylab=ylab, xlab=xlab, main=if(is.null(main)) paste("Cohort", cohorts[iplot]) else main, ...)
		else lines(this.x, this.data["median",], type="l", col=col)
		for(ipi in 1:length(pi)) {
			qs <- c((1-pi[ipi]/100)/2, (1+pi[ipi]/100)/2)
			if(! as.character(qs[1]) %in% rownames(this.data)) {
				warning("PI ", pi, " not availalble in cohort.data.")
				next
			}
			lines(this.x, this.data[as.character(qs[1]),], col=col, lty=1+ipi)
			lines(this.x, this.data[as.character(qs[2]),], col=col, lty=1+ipi)
		}
		grid()
		if(show.legend && ann)
			legend(legend.pos, legend=legend, lty=lty, bty='n', col=cols)
	}
	if(nplots > 1) par(mgp=cur.mgp, mar=cur.mar, oma=cur.oma)
}
