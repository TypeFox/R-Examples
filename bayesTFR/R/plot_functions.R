DLcurve.plot.all <- function (mcmc.list = NULL, sim.dir = NULL, 
					output.dir = file.path(getwd(), "DLcurves"), 
					output.type="png",
					burnin = NULL, verbose = FALSE,  ...) {
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	if(is.null(mcmc.list)) mcmc.list <- get.tfr.mcmc(sim.dir=sim.dir, verbose=verbose, burnin=burnin)
	mcmc.list <- get.mcmc.list(mcmc.list)
	meta <- mcmc.list[[1]]$meta
	.do.plot.all.country.loop(as.character(country.names(meta)), meta, output.dir, DLcurve.plot, output.type=output.type, 
		file.prefix='DLplot', plot.type='DL graph', verbose=verbose, mcmc.list = mcmc.list, burnin = burnin, ...)
}

stop.if.country.not.DL <- function(country.obj, meta) {
	if (!is.element(country.obj$index, meta$id_DL))
    	stop('Country ', country.obj$name, ' not estimated because no decline observed.')
}

tfr.get.dlcurves <- function(x, mcmc.list, country.code, country.index, burnin, nr.curves, predictive.distr=FALSE,
                                                               return.sigma=FALSE) {
	# if country.code is null, get world distribution
	.get.sig <- function(i, sigma, S, a, b) return(sigma + (x[i] - S)*ifelse(x[i] > S, -a, b))
	.get.sig.distr <- function(i, traces)
						return(apply(traces, 1, function(y) 
							return(pmax(.get.sig(i, y['sigma0'], y['S_sd'], y['a_sd'], y['b_sd']), mcmc$meta$sigma0.min))))
    dlc <- sigma.all <- c()
    cspec <- TRUE
    if(!is.null(country.code) && !is.element(country.index, mcmc.list[[1]]$meta$id_Tistau)) {
    	U.var <- paste("U_c", country.code, sep = "")
    	d.var <- paste("d_c", country.code, sep = "")
    	Triangle_c4.var <- paste("Triangle_c4_c", country.code, sep = "")
    	gamma.vars <- paste("gamma_", 1:3, "_c", country.code, sep = "")
    } else {
    	U.var <- "U"
    	d.var <- "d"
    	Triangle_c4.var <- "Triangle_c4"
    	gamma.vars <- paste("gamma_", 1:3, sep = "")
    	alpha.vars <- paste('alpha_',1:3, sep='')
		delta.vars <- paste('delta_',1:3, sep='')
		if(!is.null(country.code))
			Uvalue = get.observed.tfr(country.index, mcmc.list[[1]]$meta, 
										'tfr_matrix_all')[mcmc.list[[1]]$meta$tau_c[country.index]]
		else Uvalue <- mcmc.list[[1]]$meta$U.c.low.base + (mcmc.list[[1]]$meta$U.up - mcmc.list[[1]]$meta$U.c.low.base)/2
    	cspec <- FALSE
    }
    # Compute the quantiles on a sample of at least 2000.   
    nr.curves.from.mc <- if (!is.null(nr.curves)) ceiling(max(nr.curves, 2000)/length(mcmc.list))
    						else NULL
    for (mcmc in mcmc.list) {
    	th.burnin <- get.thinned.burnin(mcmc,burnin)
    	thincurves.mc <- get.thinning.index(nr.points=nr.curves.from.mc, 
            all.points=mcmc$length - th.burnin)
        if(cspec) # country specific
        	traces <- data.frame(load.tfr.parameter.traces.cs(mcmc, country.code, 
        						burnin=th.burnin, 
								thinning.index=thincurves.mc$index))
		else { #Tistau country. Use hierarchical distr.
			traces <- data.frame(load.tfr.parameter.traces(mcmc, 
        						burnin=th.burnin, 
								thinning.index=thincurves.mc$index))
			ltraces <- nrow(traces)
			for (i in 1:3){
				gamma <- rnorm(ltraces, mean = traces[,alpha.vars[i]], 
 										sd = traces[,delta.vars[i]])
				traces <- cbind(traces, gamma)
				colnames(traces)[ncol(traces)] <- gamma.vars[i]
			}
			Triangle4_tr_s <- rnorm(ltraces, mean = traces[,'Triangle4'], sd = traces[, 'delta4'])
			traces <- cbind(traces, 
						Triangle_c4=( mcmc$meta$Triangle_c4.up*exp(Triangle4_tr_s) + mcmc$meta$Triangle_c4.low)/(1+exp(Triangle4_tr_s)))
			d_tr_s <- rnorm(ltraces, mean = traces[,'chi'], sd = traces[,'psi'])
			traces <- cbind(traces,   d=(mcmc$meta$d.up*(exp(d_tr_s) + mcmc$meta$d.low)/(1+exp(d_tr_s))),
									U=rep(Uvalue, ltraces))
		}
        theta <- (traces[, U.var] - traces[, Triangle_c4.var] ) * 
            exp(traces[, gamma.vars, drop=FALSE])/apply(exp(traces[,gamma.vars, drop=FALSE]), 1, sum)
        theta <- cbind(theta, traces[, Triangle_c4.var], traces[, d.var])
        dl <- t(apply(theta, 1, DLcurve, tfr = x, p1 = mcmc$meta$dl.p1, p2 = mcmc$meta$dl.p2))
        #stop('')
        if(length(x) == 1) dl <- t(dl)
        if(predictive.distr || return.sigma) {
			wp.traces <- load.tfr.parameter.traces(mcmc, 
        						burnin=th.burnin, 
								thinning.index=thincurves.mc$index,
								par.names=c('sigma0', 'S_sd', 'a_sd', 'b_sd'))
			sigma_eps <- NULL 
			for(j in 1:length(x)) 
				sigma_eps <- cbind(sigma_eps, .get.sig.distr(j, wp.traces))
			if(predictive.distr) {
				errors <- t(apply(sigma_eps, 1, function(sig) rnorm(dim(dl)[2],0,sig)))
				dlc <- rbind(dlc, dl+errors)
			} else {
				dlc <- rbind(dlc, dl)
				sigma.all <- rbind(sigma.all, sigma_eps)
			}
        } else dlc <- rbind(dlc, dl)
    }    
    return (if(!return.sigma) dlc else list(dl=dlc, sigma=sigma.all))
}

.match.colors.with.default <- function(col, default.col) {
	ldcol <- length(default.col)
	lcol <- length(col)
	if(lcol >= ldcol) return(col)	
	col <- rep(col, ldcol)
	col[(lcol+1):ldcol] <- default.col[(lcol+1):ldcol]
	return(col)
}

DLcurve.plot <- function (mcmc.list, country, burnin = NULL, pi = 80, tfr.max = 10, 
    nr.curves = NULL, predictive.distr=FALSE, ylim = NULL, xlab = "TFR (reversed)", ylab = "TFR decrement", 
    main = NULL, show.legend=TRUE, col=c('black', 'red', "#00000020"), ...
    ) 
{	
	if(class(mcmc.list) == 'bayesTFR.prediction') {
		if(!is.null(burnin) && burnin != mcmc.list$burnin)
			warning('Prediction was generated with different burnin. Burnin set to ', mcmc.list$burnin)
		burnin <- 0 # because burnin was already cut of the traces
	}
	col <- .match.colors.with.default(col, c('black', 'red', "#00000020"))
	if(is.null(burnin)) burnin <- 0
    mcmc.list <- get.mcmc.list(mcmc.list)
    meta <- mcmc.list[[1]]$meta
    country <- get.country.object(country, meta)
    #stop.if.country.not.DL(country, meta)
    tfr_plot <- seq(0, tfr.max, 0.1)
    dlc <- tfr.get.dlcurves(tfr_plot, mcmc.list, country$code, country$index, burnin, nr.curves, 
    						predictive.distr=predictive.distr)
    miny <- min(dlc)
    maxy <- max(dlc)
    obs.data <- get.observed.tfr(country$index, meta, 'tfr_matrix_observed', 'tfr_matrix_all')[1:meta$T_end_c[country$index]]
    decr <- -diff(obs.data)
    dl.obs.idx <- if(max(meta$tau_c[country$index],1) >= meta$lambda_c[country$index]) c()
    				else seq(max(meta$tau_c[country$index],1), meta$lambda_c[country$index]-1)
    p3.obs.idx <- if(meta$lambda_c[country$index]>length(decr)) c() else seq(meta$lambda_c[country$index], length(decr))
    maxy <- max(maxy, decr, na.rm=TRUE)
    miny <- min(miny, decr, na.rm=TRUE)
    thincurves <- get.thinning.index(nr.curves, dim(dlc)[1])
    ltype <- "l"
    if (thincurves$nr.points == 0) {
        ltype <- "n"
        thincurves$index <- 1
    }
    if (is.null(main)) main <- country$name
    if (is.null(ylim)) ylim <- c(miny, maxy)
    # plot trajectories
    plot(dlc[thincurves$index[1], ] ~ tfr_plot, col = col[3], 
        type = ltype, 
        xlim = c(max(tfr_plot), min(tfr_plot)), 
        ylim = ylim, ylab = ylab, xlab = xlab, main = main, ...
        )
    if (thincurves$nr.points > 1) {
        for (i in 2:thincurves$nr.points) {
            lines(dlc[thincurves$index[i], ] ~ tfr_plot, col = col[3])
        }
    }
    # plot quantiles
    dl50 <- apply(dlc, 2, quantile, 0.5)
    lines(dl50 ~ tfr_plot, col = col[2], lwd = 2)
    lty <- 2:(length(pi) + 1)
    for (i in 1:length(pi)) {
        al <- (1 - pi[i]/100)/2
        dlpi <- apply(dlc, 2, quantile, c(al, 1 - al))
        lines(dlpi[1, ] ~ tfr_plot, col = col[2], lty = lty[i], 
            lwd = 2)
        lines(dlpi[2, ] ~ tfr_plot, col = col[2], lty = lty[i], 
            lwd = 2)
    }
    # plot observed data
    obs.legend <- list(legend=c(), pch=c(), bg=c(), col=c())
    if(length(dl.obs.idx) > 0 && any(!is.na(decr[dl.obs.idx]))) {
    	points(decr[dl.obs.idx] ~ obs.data[-length(obs.data)][dl.obs.idx], pch = 19, col=col[1])
    	obs.legend$legend <- c(obs.legend$legend, 'Phase II data')
    	obs.legend$pch <- c(obs.legend$pch, 19)
    	obs.legend$col <- c(obs.legend$col, col[1])
    	#obs.legend$bg <- c(obs.legend$bg, col[1])
    }
    endpI <- if(length(dl.obs.idx)==0) max(meta$tau_c[country$index],1) else dl.obs.idx[1]-1
    if((endpI > 1) && any(!is.na(decr[1:endpI])))  {# draw phase I points
    	points(decr[1:endpI] ~ obs.data[-length(obs.data)][1:endpI], pch=0, col=col[1])
    	obs.legend$legend <- c(obs.legend$legend, 'Phase I data')
    	obs.legend$pch <- c(obs.legend$pch, 0)
    	obs.legend$col <- c(obs.legend$col, col[1])
    	#obs.legend$bg <- c(obs.legend$bg, col[1])
    }
    if(length(p3.obs.idx)>0) {
    	points(decr[p3.obs.idx] ~ obs.data[-length(obs.data)][p3.obs.idx], pch = 2, col=col[1]) # draw phase III points
    	obs.legend$legend <- c(obs.legend$legend, 'Phase III data')
    	obs.legend$pch <- c(obs.legend$pch, 2)
    	obs.legend$col <- c(obs.legend$col, col[1])
    	#obs.legend$bg <- c(obs.legend$bg, 'grey')
    }
    if(show.legend)
    	legend("topright", legend = c("median", paste("PI", pi), obs.legend$legend), 
        	lty = c(1, lty, rep(0,length(obs.legend$pch))), bty = "n", 
        	col = c(rep(col[2], 1+length(pi)), obs.legend$col), 
        	pch=c(rep(-1,length(pi)+1), obs.legend$pch),
        	#bg=c(rep(-1,length(pi)+1), obs.legend$bg),
        	)
}

.get.trajectories.table <- function(tfr.pred, country, obs.data, pi, pred.median, cqp, half.child.variant=FALSE) {
	l <- tfr.pred$nr.projections
	obs.data <- obs.data[!is.na(obs.data)]
	x1 <- as.integer(names(obs.data))
	x2 <- seq(max(x1)+5, by=5, length=l)
	tfr <- as.matrix(obs.data, ncol=1)
	rownames(tfr) <- x1
	pred.table <- matrix(NA, ncol=2*length(pi)+1, nrow=l)
	pred.table[,1] <- pred.median[2:(l+1)]
	colnames(pred.table) <- c('median', rep(NA,ncol(pred.table)-1))
	idx <- 2
	for (i in 1:length(pi)) {
		al <- (1-pi[i]/100)/2
		if (!is.null(cqp[[i]])) {
			pred.table[,idx:(idx+1)] <- t(cqp[[i]][,2:(l+1)])
		} else{
			pred.table[,idx:(idx+1)] <- matrix(NA, nrow=l, ncol=2)
		}
		colnames(pred.table)[idx:(idx+1)] <- c(al, 1-al)
		idx <- idx+2
	}
	rownames(pred.table) <- x2
	cn <- colnames(pred.table)[2:ncol(pred.table)]
	pred.table[,2:ncol(pred.table)] <- pred.table[,cn[order(cn)]]
	colnames(pred.table)[2:ncol(pred.table)] <- cn[order(cn)]
	if(half.child.variant) {
		up.low <- get.half.child.variant(median=c(0, pred.table[,1]))
		pred.table <- cbind(pred.table, t(up.low[,2:ncol(up.low)]))
		colnames(pred.table)[(ncol(pred.table)-1):ncol(pred.table)] <- c('-0.5child', '+0.5child')
	}
	return(rbind(cbind(tfr, matrix(NA, nrow=nrow(tfr), ncol=ncol(pred.table)-1)), pred.table))
}

tfr.trajectories.table <- function(tfr.pred, country, pi=c(80, 95), half.child.variant=TRUE) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	country <- get.country.object(country, tfr.pred$mcmc.set$meta)
	obs.data <- get.data.imputed.for.country(tfr.pred, country$index)
	if(!is.null(tfr.pred$present.year.index)) obs.data <- obs.data[1:min(length(obs.data), tfr.pred$present.year.index.all)]
	pred.median <- get.median.from.prediction(tfr.pred, country$index, country$code)
	trajectories <- get.trajectories(tfr.pred, country$code)
	cqp <- list()
	for (i in 1:length(pi))
		cqp[[i]] <- get.traj.quantiles(tfr.pred, country$index, country$code, trajectories$trajectories, pi[i])
	return(.get.trajectories.table(tfr.pred, country, obs.data, pi, pred.median, cqp, half.child.variant))
}

get.typical.trajectory.index <- function(trajectories) {
	med <- apply(trajectories, 1, median, na.rm=TRUE)
	sumerrors <- apply(abs(trajectories - med), 2, sum)
	sorterrors <- order(sumerrors)
	return(sorterrors[round(length(sorterrors)/2, 0)])
}

get.trajectories <- function(tfr.pred, country, nr.traj=NULL, adjusted=TRUE, base.name='traj', typical.trajectory=FALSE) {
	traj.file <- file.path(tfr.pred$output.directory, paste(base.name, '_country', country, '.rda', sep=''))
	if (!file.exists(traj.file)) return(list(trajectories=NULL))
	load(traj.file)
	if(typical.trajectory) {
		traj.idx <- get.typical.trajectory.index(trajectories)
	} else {
		thintraj <- get.thinning.index(nr.traj, dim(trajectories)[2]) 
		if (thintraj$nr.points == 0) return(list(trajectories=NULL))
		traj.idx <- thintraj$index
	}
	if(!is.null(trajectories)) {
		if(adjusted) {
			shift <- get.tfr.shift(country, tfr.pred)
	 		if(!is.null(shift)) trajectories <- trajectories + shift
	 	}
	 	rownames(trajectories) <- get.prediction.years(tfr.pred$mcmc.set$meta, dim(trajectories)[1], 
	 													present.year.index=tfr.pred$present.year.index)
	 }
	return(list(trajectories=trajectories, index=traj.idx))
}

get.quantile.from.prediction <- function(tfr.pred, quantile, country.index, country.code=NULL, adjusted=TRUE) {
	quant.values <- tfr.pred$quantiles[country.index, as.character(quantile),]
	if (!adjusted) return(quant.values)
	shift <- get.tfr.shift(country.code, tfr.pred)
	if(!is.null(shift)) quant.values <- quant.values + shift
	return(quant.values)
}
get.median.from.prediction <- function(tfr.pred, country.index, country.code=NULL, adjusted=TRUE) {
	return(get.quantile.from.prediction(tfr.pred, quantile=0.5, country.index=country.index, 
										country.code=country.code, adjusted=adjusted))
}
	
get.traj.quantiles <- function(tfr.pred, country.index, country.code, trajectories=NULL, pi=80, adjusted=TRUE) {
	al <- (1-pi/100)/2
	quantile.values <- as.numeric(dimnames(tfr.pred$quantiles)[[2]])
	alidx<-round(quantile.values,6)==round(al,6)
	cqp <- NULL
	if (any(alidx)) { # pre-saved quantiles
		alidx2 <- round(quantile.values,6)==round(1-al,6)
		cqp <- rbind(tfr.pred$quantiles[country.index, alidx,], 
							tfr.pred$quantiles[country.index, alidx2,])
	} else { # non-standard quantiles
		reload <- FALSE
		if (is.null(trajectories)) {
			if(tfr.pred$nr.traj > 0) reload <- TRUE
		} else { 
			if (dim(trajectories)[2] < 2000 && tfr.pred$nr.traj > dim(trajectories)[2]) reload <- TRUE
		}
		if(reload) {
			#load 2000 trajectories maximum for computing quantiles
			traj.reload <- get.trajectories(tfr.pred, tfr.pred$mcmc.set$meta$regions$country_code[country.index], 2000)
			trajectories <- traj.reload$trajectories
		}
		if (!is.null(trajectories)) {
			cqp <- apply(trajectories, 1, 
						quantile, c(al, 1-al), na.rm = TRUE)
		}
	}
	if(!adjusted) return(cqp)
	shift <- get.tfr.shift(country.code, tfr.pred)
	if(!is.null(shift)) cqp <- cqp + matrix(shift, nrow=nrow(cqp), ncol=ncol(cqp), byrow=TRUE)
	return(cqp)
}
	
tfr.trajectories.plot.all <- function(tfr.pred, 
									output.dir=file.path(getwd(), 'TFRtrajectories'),
									output.type="png", main=NULL, verbose=FALSE, ...) {
	# plots TFR trajectories for all countries
	.do.plot.all(tfr.pred$mcmc.set$meta, output.dir, tfr.trajectories.plot, output.type=output.type, 
						verbose=verbose, tfr.pred=tfr.pred, ...)
}

.do.plot.all.country.loop <- function(country.names, meta, output.dir, func, output.type="png", 
						file.prefix='TFRplot', plot.type='TFR graph', country.table=NULL,
						main=NULL, verbose=FALSE, ...) {					
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	main.arg <- main
	for (country in country.names) {
		country.obj <- if(!is.null(meta)) get.country.object(country, meta)
						else get.country.object(country, country.table=country.table)
		if(verbose)
			cat('Creating', plot.type, 'for', country, '(', country.obj$code, ')\n')
		if(!is.null(main) && grepl('XXX', main, fixed=TRUE))
			main.arg <- gsub('XXX', as.character(country.obj$name), main, fixed=TRUE)
		do.call(output.type, list(file.path(output.dir, 
										paste(file.prefix,'_c', country.obj$code, '.', postfix, sep=''))))
		do.call(func, list(country=country.obj$code, main=main.arg, ...))
		dev.off()
	}
	if(verbose)
		cat('\nPlots stored into', output.dir, '\n')	
}

.do.plot.all <- function(meta, ...) {
	# processes plotting function func for all countries
	.do.plot.all.country.loop(country.names(meta), meta, ...)				
}


get.half.child.variant <- function(median, increment=c(0, 0.25, 0.4, 0.5)) {
	l <- length(median)
	lincr <- length(increment)
	upper <- lower <- c()
	for (i in 1:l) {
		upper <- c(upper, median[i]+increment[min(i,lincr)])
		lower <- c(lower, median[i]-increment[min(i,lincr)])
	}
	return(rbind(lower, upper))	
}

tfr.trajectories.plot <- function(tfr.pred, country, pi=c(80, 95), 
								  half.child.variant=TRUE, nr.traj=NULL,
								  adjusted.only = TRUE, typical.trajectory=FALSE,
								  mark.estimation.points=FALSE,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='Year', ylab='TFR', main=NULL, lwd=c(2,2,2,2,2,1), 
								  col=c('black', 'green', 'red', 'red', 'blue', '#00000020'),
								  show.legend=TRUE, add=FALSE, ...
								  ) {
	# lwd/col is a vector of 6 line widths/colors for: 
	#	1. observed data, 2. imputed missing data, 3. median, 4. quantiles, 5. +- 0.5 child, 6. trajectories
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	if(length(lwd) < 6) {
		lwd <- rep(lwd, 6)
		lwd[6] <- 1
	}
	col <- .match.colors.with.default(col, c('black', 'green', 'red', 'red', 'blue', '#00000020'))
	country <- get.country.object(country, tfr.pred$mcmc.set$meta)
	tfr_observed <- get.observed.tfr(country$index, tfr.pred$mcmc.set$meta, 'tfr_matrix_observed', 'tfr_matrix_all')
	T_end_c <- tfr.pred$mcmc.set$meta$T_end_c
	if(!is.null(tfr.pred$present.year.index.all)) T_end_c <- pmin(T_end_c, tfr.pred$present.year.index.all)
	tfr_matrix_reconstructed <- get.tfr.reconstructed(tfr.pred$tfr_matrix_reconstructed, tfr.pred$mcmc.set$meta)
	suppl.T <- length(tfr_observed) - tfr.pred$mcmc.set$meta$T_end
	y1.part1 <- tfr_observed[1:T_end_c[country$index]]
	y1.is.not.na <- which(!is.na(y1.part1))
	y1.part1 <- y1.part1[y1.is.not.na]
	lpart1 <- length(y1.part1)
	y1.part2 <- NULL
	lpart2 <- min(tfr.pred$mcmc.set$meta$T_end, tfr.pred$present.year.index) - T_end_c[country$index] + suppl.T
	if (lpart2 > 0) {
		p2idx <- (T_end_c[country$index]+1-suppl.T):nrow(tfr_matrix_reconstructed)
		y1.part2 <- tfr_matrix_reconstructed[p2idx,country$index]
		names(y1.part2) <- rownames(tfr_matrix_reconstructed)[p2idx]
	}
	x1 <- as.integer(c(names(y1.part1), names(y1.part2)))
	x2 <- as.numeric(dimnames(tfr.pred$quantiles)[[3]])
	trajectories <- get.trajectories(tfr.pred, country$code, nr.traj, typical.trajectory=typical.trajectory)
	# plot historical data: observed
	if (!add) {
		if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
		if(is.null(ylim)) ylim <- c(0, max(trajectories$trajectories, y1.part1, y1.part2, na.rm=TRUE))
		if(is.null(main)) main <- country$name
		plot(xlim, ylim, type='n', xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
					panel.first = grid())
	}
	points.x <- x1[1:lpart1]
	points.y <- y1.part1
	if (mark.estimation.points) {
		tfr.est <- get.observed.tfr(country$index, tfr.pred$mcmc.set$meta, 'tfr_matrix')[1:T_end_c[country$index]][y1.is.not.na]
		end.na <- which(!is.na(tfr.est))
		end.na <- if(length(end.na)==0) length(tfr.est) else end.na[1]
		if(end.na > 1) {
			na.idx <- 1:end.na
			points(points.x[na.idx], points.y[na.idx], type=type, lwd=lwd[1], col=rgb(t(col2rgb(col[1])/255), alpha=0.1), ...)
			points.x <- points.x[-na.idx[-end.na]]
			points.y <- points.y[-na.idx[-end.na]]
		}
	}
	points(points.x, points.y, type=type, lwd=lwd[1], col=col[1], ...)
	if(lpart2 > 0) { # imputed values
		lines(x1[(lpart1+1): length(x1)], y1.part2, pch=2, type='b', col=col[2], lwd=lwd[2])
		lines(x1[lpart1:(lpart1+1)], c(y1.part1[lpart1], y1.part2[1]), col=col[2], lwd=lwd[2]) # connection between the two parts
	}
	
	# plot trajectories
	if(!is.null(trajectories$trajectories)) { 
		for (i in 1:length(trajectories$index)) {
			lines(x2, trajectories$trajectories[,trajectories$index[i]], type='l', col=col[6], lwd=lwd[6])
		}
	}
	# plot median
	tfr.median <- get.median.from.prediction(tfr.pred, country$index, country$code)
	lines(x2, tfr.median, type='l', col=col[3], lwd=lwd[3]) 
	# plot given CIs
	lty <- 2:(length(pi)+1)
	for (i in 1:length(pi)) {
		cqp <- get.traj.quantiles(tfr.pred, country$index, country$code, trajectories$trajectories, pi[i])
		if (!is.null(cqp)) {
			lines(x2, cqp[1,], type='l', col=col[4], lty=lty[i], lwd=lwd[4])
			lines(x2, cqp[2,], type='l', col=col[4], lty=lty[i], lwd=lwd[4])
		}
	}
	legend <- c()
	cols <- c()
	lwds <- c()
	lty <- c(1, lty)
	if(!adjusted.only) { # plot unadjusted median
		bhm.median <- get.median.from.prediction(tfr.pred, country$index, country$code, adjusted=FALSE)
		lines(x2, bhm.median, type='l', col=col[3], lwd=lwd[3], lty=max(lty)+1)
		legend <- c(legend, 'BHM median')
		cols <- c(cols, col[3])
		lwds <- c(lwds, lwd[3])
		lty <- c(max(lty)+1, lty)
	}
	median.legend <- if(adjusted.only) 'median' else 'adj. median'
	legend <- c(legend, median.legend, paste(pi, '% PI', sep=''))
	cols <- c(cols, col[3], rep(col[4], length(pi)))
	lwds <- c(lwds, lwd[3], rep(lwd[4], length(pi)))
	if (half.child.variant) {
		lty <- c(lty, max(lty)+1)
		llty <- length(lty)
		up.low <- get.half.child.variant(median=tfr.median)
		lines(x2, up.low[1,], type='l', col=col[5], lty=lty[llty], lwd=lwd[5])
		lines(x2, up.low[2,], type='l', col=col[5], lty=lty[llty], lwd=lwd[5])
		legend <- c(legend, '+/- 0.5 child')
		cols <- c(cols, col[5])
		lwds <- c(lwds, lwd[5])
	}
	if(show.legend) {
		legend <- c(legend, 'observed TFR')
		cols <- c(cols, col[1])
		lty <- c(lty, 1)
		pch <- c(rep(-1, length(legend)-1), 1)
		lwds <- c(lwds, lwd[1])
		if(lpart2 > 0) {
			legend <- c(legend, 'imputed TFR')
			cols <- c(cols, col[2])
			lty <- c(lty, 1)
			pch <- c(pch, 2)
			lwds <- c(lwds, lwd[2])
		}
		legend('bottomleft', legend=legend, lty=lty, bty='n', col=cols, pch=pch, lwd=lwds)
	}
}

extract.plot.args <- function(...) {
	# split '...' into plot arguments and the rest
	all.plot.args <- names(formals(plot.default))
	args <- list(...)
	which.plot.args <- pmatch(names(args), all.plot.args)
	is.fun.arg <- is.na(which.plot.args)
	return(list(plot.args=args[!is.fun.arg], other.args=args[is.fun.arg]))
}
	
do.plot.tfr.partraces <- function(mcmc.list, func, par.names, main.postfix='', chain.ids=NULL, 
									nr.points=NULL, dev.ncol=5, ...) {
	mcmc.list <- get.mcmc.list(mcmc.list)
	if (is.null(chain.ids)) {
		nr.chains <- length(mcmc.list)
		chain.ids <- 1:nr.chains
	} else {
		nr.chains <- length(chain.ids)
		}
	pars <- list()
	iter <- rep(NA, nr.chains)
	mclen <- rep(0, nr.chains)
	
	# split '...' into function arguments and plot arguments
	split.args <- extract.plot.args(...)
	fun.args <- split.args$other.args
	plot.args <- split.args$plot.args
	thin <- fun.args$thin
	fun.args$thin <- NULL
	if(is.null(fun.args$burnin)) fun.args$burnin <- 0
	orig.burnin <- fun.args$burnin
	i <- 1
	for(chain in chain.ids) {
		mcmc <- mcmc.list[[chain]]
		if (!is.null(thin) || mcmc$thin > 1) {
			consolidated.burn.thin <- burn.and.thin(mcmc, orig.burnin, 
										if (is.null(thin)) mcmc$thin else thin)
			fun.args$burnin <- consolidated.burn.thin$burnin
			if (!is.null(consolidated.burn.thin$index)) fun.args$thinning.index <- consolidated.burn.thin$index
			else {
				if(!is.null(thin)) {
					thin <- max(thin, mcmc$thin)
					fun.args$thinning.index <- unique(round(seq(1,mcmc$length-consolidated.burn.thin$burnin, 
												by=thin/mcmc$thin)))
				} else  {
						fun.args$thinning.index <- seq(1,mcmc$length-consolidated.burn.thin$burnin)
						mclen[i] <- length(fun.args$thinning.index)
					}
			}
			
		}
		pars[[i]] <- eval(do.call(func, c(list(mcmc, par.names=par.names), fun.args)))
		pars[[i]] <- filter.traces(pars[[i]], par.names)
		iter[i] <- mcmc$finished.iter
		if (i==1) {
			par.names.l <- length(colnames(pars[[1]]))
			maxy <- rep(NA, par.names.l)
			miny <- rep(NA, par.names.l)
		}
		ipara<-1
		for (para in colnames(pars[[i]])) {
			maxy[ipara] <- max(maxy[ipara], pars[[i]][,para], na.rm=TRUE)
			miny[ipara] <- min(miny[ipara], pars[[i]][,para], na.rm=TRUE)
			ipara <- ipara+1
		}
		i <- i+1
	}

	if (par.names.l < dev.ncol) {
		ncols <- par.names.l
		nrows <- 1
	} else {
		ncols <- dev.ncol
		nrows <- ceiling(par.names.l/dev.ncol)
	}
	par(mfrow=c(nrows,ncols))
	col <- 1:nr.chains
	maxx<-max(iter)
	if(is.null(plot.args$xlim)) plot.args$xlim <- c(1+orig.burnin,maxx)
	if(is.null(plot.args$xlab)) plot.args$xlab <- 'iterations'
	if(is.null(plot.args$ylab)) plot.args$ylab <- ''
	ylim <- plot.args$ylim
	ipara <- 1
	for (para in colnames(pars[[1]])) {
		#mx <- length(pars[[1]][,para])
		if (is.null(thin)) {
			maxmclen <- max(mclen)
			xindex <- if(maxmclen > 0) seq(1+orig.burnin, maxx, length=maxmclen)
						else (1+orig.burnin):maxx 
		} else xindex <- seq(1+orig.burnin, maxx, by=thin)
		thinpoints <- get.thinning.index(nr.points, length(xindex))
		if (thinpoints$nr.points > 0) {
			plot.args$ylim <- if(is.null(ylim)) c(miny[ipara], maxy[ipara]) else ylim
			do.call('plot', c(list(xindex[thinpoints$index], 
						pars[[1]][thinpoints$index, para], 
						main=paste(para, main.postfix),
						col=col[1], type='l'), plot.args))
			if (nr.chains > 1) {
				for (i in 2:nr.chains) {
					lines(xindex[thinpoints$index], 
						pars[[i]][thinpoints$index,para], col=col[i], type='l')
				}
			}
		}
		ipara <- ipara+1
	}
	#stop('')
}

tfr.partraces.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr.parameter.names(trans=TRUE), 
									nr.points=NULL, dev.ncol=5, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	do.plot.tfr.partraces(mcmc.list, 'load.tfr.parameter.traces', chain.ids=chain.ids, 
							nr.points=nr.points, par.names=par.names, dev.ncol=dev.ncol, ...)
}

tfr.partraces.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'),
									chain.ids=NULL, par.names=tfr.parameter.names.cs(trans=TRUE),
									nr.points=NULL, dev.ncol=3, low.memory=TRUE, ...) {

	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	mcmc.list <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.list[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	stop.if.country.not.DL(country.obj, mcmc.list[[1]]$meta)
	do.plot.tfr.partraces(mcmc.list, 'load.tfr.parameter.traces.cs', 
		main.postfix=paste('(',country.obj$name,')', sep=''), chain.ids=chain.ids, nr.points=nr.points, 
		country=country.obj$code, par.names=par.names, dev.ncol=dev.ncol, ...)
}

tfr3.partraces.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr3.parameter.names(), 
									nr.points=NULL, dev.ncol=3, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr3.mcmc(sim.dir, low.memory=low.memory)
	else if(class(mcmc.list)=='bayesTFR.prediction')
			stop('Function not available for bayesTFR.prediction objects.')
	tfr.partraces.plot(mcmc.list, sim.dir=NULL, chain.ids=chain.ids, par.names=par.names, 
						nr.points=nr.points, dev.ncol=dev.ncol, ...)
}

tfr3.partraces.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'),
									chain.ids=NULL, par.names=tfr3.parameter.names.cs(), 
									nr.points=NULL, dev.ncol=2, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr3.mcmc(sim.dir, low.memory=low.memory)
	else if(class(mcmc.list)=='bayesTFR.prediction')
			stop('Function not available for bayesTFR.prediction objects.')
	mcmc.list <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.list[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	do.plot.tfr.partraces(mcmc.list, 'load.tfr.parameter.traces.cs', 
		main.postfix=paste('(',country.obj$name,')', sep=''), chain.ids=chain.ids, nr.points=nr.points, 
		country=country.obj$code, par.names=par.names, dev.ncol=dev.ncol, ...)
}
		
do.plot.tfr.pardensity <- function(mcmc.list, func, par.names, par.names.ext, main.postfix='', 
								func.args=NULL, chain.ids=NULL, burnin=NULL, dev.ncol=5, ...) {
	if(class(mcmc.list) == 'bayesTFR.prediction') {
		if(!is.null(burnin) && burnin != mcmc.list$burnin)
			warning('Prediction was generated with different burnin. Burnin set to ', mcmc.list$burnin,
					'.\n Use a bayesTFR.mcmc.set object as the first argument, if the original traces should be used.')
		burnin <- 0 # because burnin was already cut of the traces
		if (!is.null(chain.ids) && max(chain.ids) > 1) {
			warning('Thinned traces from all chains used for plotting density.\n For selecting individual chains, use a bayesTFR.mcmc.set object as the first argument.')
			chain.ids <- NULL
		}
	}
	if(is.null(burnin)) burnin <- 0

	mcmc.list <- get.mcmc.list(mcmc.list)
	if (!is.null(chain.ids)) mcmc.list <- mcmc.list[chain.ids]
	par.names.l <- length(par.names.ext)
	if (par.names.l < dev.ncol) {
		ncols <- par.names.l
		nrows <- 1
	} else {
		ncols <- dev.ncol
		nrows <- ceiling(par.names.l/dev.ncol)
	}
	args <- extract.plot.args(...)
	par(mfrow=c(nrows,ncols))
	for (para in par.names) {
		values <- eval(do.call(func, c(list(mcmc.list, par.names=para, burnin=burnin), func.args)))
		values <-  filter.traces(values, par.names)
		for (par.name in colnames(values)) {
			dens <- do.call('density', c(list(values[,par.name]), args$other.args))
			do.call('plot', c(list(dens, main=paste(par.name, main.postfix)), args$plot.args))
		}
	}
}

tfr.pardensity.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr.parameter.names(trans=TRUE), 
									burnin=NULL, dev.ncol=5, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	par.names.ext <- get.full.par.names(par.names, tfr.parameter.names.extended())
	if(length(par.names.ext) <= 0)
		stop('Parameter names are not valid parameters.\nUse function tfr.parameter.names(...) or valid parameter names.')
	do.plot.tfr.pardensity(mcmc.list, 'get.tfr.parameter.traces', chain.ids=chain.ids, par.names=par.names,
							par.names.ext=par.names.ext,
							burnin=burnin, dev.ncol=dev.ncol, ...)
}

tfr.pardensity.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr.parameter.names.cs(trans=TRUE), 
									burnin=NULL, dev.ncol=3, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	mcmc.l <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.l[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	stop.if.country.not.DL(country.obj, mcmc.l[[1]]$meta)
	par.names.ext <- get.full.par.names.cs(par.names, 
											tfr.parameter.names.cs.extended(country.obj$code))
	if(length(par.names.ext) <= 0)
		stop('Parameter names are not valid country-specific parameters.\nUse function tfr.parameter.names.cs(...) or valid parameter names.')
	do.plot.tfr.pardensity(mcmc.list, 'get.tfr.parameter.traces.cs', chain.ids=chain.ids, par.names=par.names,
							par.names.ext=par.names.ext,
							main.postfix=paste('(',country.obj$name,')', sep=''),
							func.args=list(country.obj=country.obj),
							burnin=burnin, dev.ncol=dev.ncol, ...)
}

tfr3.pardensity.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr3.parameter.names(), 
									burnin=NULL, dev.ncol=3, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr3.mcmc(sim.dir, low.memory=low.memory)
	else if(class(mcmc.list)=='bayesTFR.prediction')
			stop('Function not available for bayesTFR.prediction objects.')
	do.plot.tfr.pardensity(mcmc.list, 'get.tfr.parameter.traces', chain.ids=chain.ids, par.names=par.names,
							par.names.ext=par.names,
							burnin=burnin, dev.ncol=dev.ncol, ...)
}

tfr3.pardensity.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr3.parameter.names.cs(), 
									burnin=NULL, dev.ncol=2, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr3.mcmc(sim.dir, low.memory=low.memory)
	else if(class(mcmc.list)=='bayesTFR.prediction')
			stop('Function not available for bayesTFR.prediction objects.')
	mcmc.l <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.l[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	par.names.ext <- get.full.par.names.cs(par.names, paste(par.names, '_c', country.obj$code, sep=''))
	if(length(par.names.ext) <= 0)
		stop('Parameter names are not valid country-specific parameters.\nUse function tfr3.parameter.names.cs(...) or valid parameter names.')
	do.plot.tfr.pardensity(mcmc.list, 'get.tfr.parameter.traces.cs', chain.ids=chain.ids, par.names=par.names,
							par.names.ext=par.names.ext,
							main.postfix=paste('(',country.obj$name,')', sep=''),
							func.args=list(country.obj=country.obj),
							burnin=burnin, dev.ncol=dev.ncol, ...)
}

".get.gamma.pars" <- function(pred, ...) UseMethod (".get.gamma.pars")
.get.gamma.pars.bayesTFR.prediction <- function(pred, ...) {
	# estimated by
	# library(MASS)
	# data <- pred$tfr_matrix_reconstructed[12,]
	# gd <- fitdistr(data-min(data)+0.05, densfun='gamma')
	# min(data) is 0.95
	return(list(gamma.pars=list(shape=1.87, rate=0.94), gamma.shift=0.95-0.05, min.value=0.7,
					max.value=NULL))
}

get.tfr.map.parameters <- function(pred, tfr.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, ...) {
	map.pars <- list(pred=pred, quantile=quantile, ...)
	if (same.scale) {
		gp <- .get.gamma.pars(pred)
		data <- pred$quantiles[,as.character(quantile),1]
		q <- if(is.null(tfr.range)) c(min(pmax(data,gp$gamma.shift)), max(data))-gp$gamma.shift
			 else c(max(gp$gamma.shift, tfr.range[1]-gp$gamma.shift), max(gp$gamma.shift, tfr.range[2]-gp$gamma.shift))
		min.max.q <- pgamma(q, shape=gp$gamma.pars[['shape']], rate=gp$gamma.pars[['rate']])
		quantiles <- seq(min.max.q[1], min.max.q[2], length=nr.cats)
		quant.values <- c(gp$min.value, 
				qgamma(quantiles, shape=gp$gamma.pars[['shape']], rate=gp$gamma.pars[['rate']])+gp$gamma.shift)
		if(!is.null(gp$max.value) && gp$max.value > max(quant.values)) quant.values <- c(quant.values, gp$max.value)

		if(!is.null(tfr.range)) {
			if(tfr.range[1] < quant.values[1])
				quant.values <- c(tfr.range[1], quant.values)
			if(tfr.range[1] > quant.values[1])
				quant.values <- quant.values[quant.values >= tfr.range[1]]
			last <- length(quant.values)
			if(tfr.range[2] > quant.values[last])
				quant.values <- c(quant.values, tfr.range[2])
			if(tfr.range[2] < quant.values[last])
				quant.values <- quant.values[quant.values <= tfr.range[2]]
		}
		col <- rainbow(500, start=0, end=0.67)[seq(500, 1, length=length(quant.values)-1)]
		map.pars$catMethod <- quant.values
	} else {
		col <- rainbow(500, start=0, end=0.67)[seq(500, 1, length=nr.cats)]
		map.pars$numCats <- nr.cats
	}
	map.pars$colourPalette <- col
	return(map.pars)		
}

bdem.map.all <- function(pred, output.dir, type='tfr', output.type='png', range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, file.prefix='TFRwrldmap_', ...) {
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.years <- get.all.prediction.years(pred)
	output.args <- list()
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	filename.arg <- 'filename'
	if(output.type=='postscript' || output.type=='pdf') {filename.arg<-'file'}
	else{output.args[['width']] <- 1000}
	
	map.pars <- get.tfr.map.parameters(pred, tfr.range=range, nr.cats=nr.cats, same.scale=same.scale,
					quantile=quantile, ...)
	
	for (year in all.years) {
		output.args[[filename.arg]] <- file.path(output.dir, 
										paste(file.prefix, year, '.', postfix, sep=''))
		do.call(paste(type, '.map', sep=''), c(list(year=year, device=output.type, 
							device.args=output.args), map.pars))
		dev.off()
	}
	cat('Maps written into', output.dir, '\n')
}

tfr.map.all <- function(pred, output.dir, output.type='png', tfr.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, file.prefix='TFRwrldmap_', ...) {
	bdem.map.all(pred=pred, output.dir=output.dir, type='tfr', output.type=output.type, range=tfr.range,
						nr.cats=nr.cats, same.scale=same.scale, quantile=quantile, file.prefix=file.prefix, ...)
}

".map.main.default" <- function(pred, ...) UseMethod (".map.main.default")

.map.main.default.bayesTFR.prediction <- function(pred, ...) return('TFR: quantile')
	
par.names.for.worldmap <- function(pred, ...) UseMethod ("par.names.for.worldmap")

par.names.for.worldmap.bayesTFR.prediction <- function(pred, ...) {
	return(c('lambda', tfr.parameter.names.cs.extended()))
}

"get.data.for.worldmap" <- function(pred, ...) UseMethod ("get.data.for.worldmap")

get.data.for.worldmap.bayesTFR.prediction <- function(pred, quantile=0.5, year=NULL, par.name=NULL, 
									adjusted=FALSE, projection.index=1, pi=NULL, ...) {
	meta <- pred$mcmc.set$meta
	quantiles <- quantile
	if (!is.null(pi)) {
		qlower <- (1-pi/100)/2
		quantiles <- c(quantile, qlower, 1-qlower)
	}
	if(!is.null(par.name)) { # data are parameter values
		if (!is.element(par.name, par.names.for.worldmap(pred)))
				stop('Illegal par.name. Allowed values:', 
						paste(par.names.for.worldmap(pred), collapse=', '))
		data <- c()
		if (par.name == 'lambda') {
			tfr <- get.data.imputed(pred)
			tfr.years <- get.estimation.years(meta)
			all.years <- c(tfr.years, get.prediction.years(meta, pred$nr.projections+1, pred$present.year.index)[-1])
			nr.data <- pred$nr.projections+dim(tfr)[1]
			for (country in 1:get.nr.countries(meta)) {
				country.obj <- get.country.object(country, meta, index=TRUE)
				tfr.and.pred.median <- c(tfr[,country], 
						get.quantile.from.prediction(pred, quantile, country.obj$index, country.obj$code, 
												adjusted=adjusted)[-1])
				lambda <- all.years[find.lambda.for.one.country(tfr.and.pred.median, nr.data)]
				data <- c(data, lambda)
			}
			codes <- meta$regions$country_code
		} else {
			con <- textConnection("sout", "w", local=TRUE) # redirect output (to get rid of coda's messages)
			for (country in get.countries.index(meta)) {
				country.obj <- get.country.object(country, meta, index=TRUE)
				sink(con, type='message')
				s <- summary(coda.list.mcmc(pred$mcmc.set, country=country.obj$code, 
						par.names=NULL, par.names.cs=par.name, thin=1, burnin=0), quantiles = quantiles)
				sink(type='message')
				data <- rbind(data, s$quantiles)
			}
			close(con)
			codes <- meta$regions$country_code[get.countries.index(meta)]
		}
		projection.index <- 1
		projection <- TRUE
		period <- paste('Parameter', par.name, 'for')
	} else { # data are TFRs
		projection <- TRUE
		if(!is.null(year)) {
			ind.proj <- get.predORest.year.index(pred, year)
			projection.index <- ind.proj['index']
			projection <- ind.proj['is.projection']
			if(is.null(projection.index)) 
				if(is.null(projection.index)) stop('Projection year ', year, ' not found.')
		}
		if(projection) {
			if(!all(is.element(as.character(quantiles), dimnames(pred$quantiles)[[2]])))
				stop('Some of the quantiles ', paste(quantiles, collapse=', '), ' not found.\nAvailable: ', 
							paste(dimnames(pred$quantiles)[[2]], collapse=', '), 
					 '\nCheck arguments "quantile" and "pi".')
			data <- pred$quantiles[, as.character(quantiles), projection.index]
			if(adjusted) data <- data + get.tfr.shift.all(pred, projection.index)
			period <- get.prediction.periods(meta, projection.index, 
								present.year.index=pred$present.year.index)[projection.index]
		} else {
			data <- get.data.imputed(pred)[projection.index, ]
			period <- get.tfr.periods(meta)[projection.index]
		}
		codes <- meta$regions$country_code
	}
	if(adjusted) period <- paste(period, 'adjusted')
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
	return(list(period=period, data=res, country.codes=codes, lower=low, upper=up))
}

tfr.map <- function(pred, quantile=0.5, year=NULL, par.name=NULL, adjusted=FALSE, 
					projection.index=1,  device='dev.new', main=NULL, 
					resolution=c("coarse","low","less islands","li","high"),
					device.args=NULL, data.args=NULL, ...
				) {
	resolution <- match.arg(resolution)
	#if(resolution=='high') require(rworldxtra)
	data.period <- do.call(get.data.for.worldmap, c(list(pred, quantile, year=year, 
									par.name=par.name, adjusted=adjusted, projection.index=projection.index), data.args))
	#data.period.base <- do.call(get.data.for.worldmap, c(list(pred, quantile, year=2013, 
	#								par.name=par.name, adjusted=adjusted, projection.index=projection.index), data.args))
	#data <- (data.period$data - data.period.base$data)/1000
	data <- data.period$data
	period <- data.period$period
	tfr <- data.frame(cbind(un=data.period$country.codes, tfr=data))
	map <- rworldmap::getMap(resolution=resolution)
	#first get countries excluding Antarctica which crashes spTransform (says the help page for joinCountryData2Map)
	sPDF <- map[-which(map$ADMIN=='Antarctica'), ]	
	if(requireNamespace("rgdal", quietly=TRUE)) {
		#transform map to the Robinson projection
		sPDF <- sp::spTransform(sPDF, CRSobj=sp::CRS("+proj=robin +ellps=WGS84"))
	}
	## recode missing UN codes and UN member states
	sPDF$UN <- sPDF$ISO_N3
	## N. Cyprus -> assign to Cyprus
	sPDF$UN[sPDF$ISO3=="CYN"] <- 196
	## Kosovo -> assign to Serbia
	sPDF$UN[sPDF$ISO3=="KOS"] <- 688
	## W. Sahara -> no UN numerical code assigned in Natural Earth map
	sPDF$UN[sPDF$ISO3=="SAH"] <- 732
	## Somaliland -> assign to Somalia
	sPDF$UN[sPDF$ISO3=="SOL"] <- 706

	#mtfr <- joinCountryData2Map(tfr, joinCode='UN', nameJoinColumn='un')
	# join sPDF with tfr
	mtfr <- rep(NA, length(sPDF$UN))
	valididx <- which(is.element(sPDF$UN, tfr$un))
	mtfr[valididx] <- tfr$tfr[sapply(sPDF$UN[valididx], function(x,y) which(y==x),  tfr$un)]
	sPDF$tfr <- mtfr
	if(is.null(main)) {
		main <- paste(period, .map.main.default(pred, data.period), quantile)
	}
	if (device != 'dev.cur')
		do.call(rworldmap::mapDevice, c(list(device=device), device.args))
	mapParams<-rworldmap::mapCountryData(sPDF, nameColumnToPlot='tfr', addLegend=FALSE, mapTitle=main, ...
	)
	do.call(rworldmap::addMapLegend, c(mapParams, legendWidth=0.5, legendMar=2, legendLabels='all'))
	#do.call(addMapLegend, c(mapParams, legendWidth=0.5, legendMar=2, legendLabels='all', sigFigs=2, legendShrink=0.8, tcl=-0.3, digits=1))
}

tfr.map.gvis <- function(pred, year=NULL, quantile=0.5, pi=80, par.name=NULL, 
							  adjusted=FALSE, ...)
	bdem.map.gvis(pred, year=year, 
						quantile=quantile, pi=pi, par.name=par.name, adjusted=adjusted, ...)


"bdem.map.gvis" <- function(pred, ...) UseMethod ("bdem.map.gvis")

bdem.map.gvis.bayesTFR.prediction <- function(pred, year=NULL, quantile=0.5, pi=80, 
										par.name=NULL, html.file=NULL, adjusted=FALSE, ...) {
	.do.gvis.bdem.map('TFR', 'BHM for Total Fertility Rate<br>', pred, year=year, 
						quantile=quantile, pi=pi, par.name=par.name, adjusted=adjusted, ...)
}

.do.gvis.bdem.map <- function(what, title, pred, year=NULL, quantile=0.5, pi=80, 
									par.name=NULL, adjusted=FALSE, ...) {
	data.period <- get.data.for.worldmap(pred, quantile, year=year, 
									par.name=par.name, projection.index=1, adjusted=adjusted, pi=pi, ...)
	mapdata <- data.period$data
	period <- data.period$period
	lower <- data.period$lower
	upper <- data.period$upper
 	un <- data.period$country.codes
 	countries.table <- get.countries.table(pred)
 	if(!is.null(par.name)) what <- par.name
 	e <- new.env()
	data('iso3166', envir=e)
	unmatch <- match(un, e$iso3166$uncode)
	unidx <- which(!is.na(unmatch))	
	ct.idx <- sapply(un[unidx], function(x, y) which(y==x), countries.table$code)
	country.names <- countries.table$name[ct.idx]
	#remove problematic characters
	country.names <- iconv(country.names, "latin1", "ASCII", "?") 
	data <- data.frame(un=un[unidx], name=country.names, 
						iso=e$iso3166$charcode[unmatch][unidx])	
	data[[what]] <- mapdata[unidx]
	if(!is.null(lower)) { # confidence intervals defined
		lower.name <- paste('lower_', pi, sep='')
		upper.name <- paste('upper_', pi, sep='')
		data[[lower.name]] <- round(lower[unidx], 2)
		data[[upper.name]] <- round(upper[unidx], 2)
		data$pi <- paste(e$iso3166$charcode[unmatch][unidx], ': ', pi, '% PI (', data[[lower.name]], ', ', 
				data[[upper.name]], ')', sep='')
		hovervar <- 'pi'
	} else { # no confidence intervals
		lower.name <- 'lower'
		upper.name <- 'upper'
		data[[lower.name]] <- data[[upper.name]] <- rep(NA, length(unidx))
		hovervar <- ''
	}
	col <- c('0x0000CC', '0x00CCFF', '0x33FF66', '0xFFFF66', '0xFF9900', '0xFF3300')
	geo <- googleVis::gvisGeoMap(data, locationvar="iso", numvar=what, hovervar=hovervar, 
				options=list(height=500, width=900, dataMode='regions',
				colors=paste('[', paste(col, collapse=', '), ']')))

    #geo$html$caption <- paste(title, 'in', period,'<br>\n')
    geo$html$caption <- paste(title, period, .map.main.default(pred, data.period), quantile)
    bdem.data <- data[,c('un', 'iso', 'name', what, lower.name, upper.name)]
	gvis.table <- googleVis::gvisTable(bdem.data, 
							options=list(width=600, height=600, page='disable', pageSize=198))
	page <- list(type="MapTable", 
 			 	 chartid=format(Sys.time(), "BdemMap-%Y-%m-%d-%H-%M-%S"), 
 			 		html=list(Header=geo$html$header,
 					Chart1=geo$html$chart,
 					Caption1=geo$html$caption,
 					Chart2=gvis.table$html$chart,
 					Caption2=gvis.table$html$caption,               
 					Footer=gvis.table$html$footer)
             )
	class(page) <- list("gvis", class(page))
	plot(page)
	invisible(page)
}