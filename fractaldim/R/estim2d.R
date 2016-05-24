has.direction <- function(all.directions, direction) {
	return(length(grep(direction, all.directions, fixed=TRUE)) > 0)
}

do.fd.estimate2d.simple <- function(data, p=2, method.name, method1d, direction='hv',
						 	plot.loglog=FALSE, nlags = "auto", 
							plot.allpoints=FALSE, legend.type='s', ...,
                          	debuglevel=0) {
# Find, for every row and every column, an estimate of D, by using the respective 1d estimator. 
# Then get the median over these estimates. 
	repl.row <- if (has.direction(direction, 'h')) dim(data)[1] else 0
	repl.col <- if (has.direction(direction, 'v')) dim(data)[2] else 0

	fd <- matrix(NA, nrow=2, ncol=repl.row+repl.col, 
					dimnames=list(c('fd', 'scale'), c()))
	args <- list(plot.loglog=FALSE, nlags=nlags,
				 plot.allpoints=FALSE, ..., 
				 debuglevel=debuglevel)
				 
	if (repl.row > 0) {
		for (i in 1:repl.row) {
			if (debuglevel>=5) 
				cat("\nRow:", i, ", method: transect\n--------------------")
			res <- do.call(paste('fd.estim.', method1d, sep=''), c(list(data[i,], p.index=p), args))
			fd[,i] <- c(res$fd, res$scale)
		}
	}
	if (repl.col > 0) {
		for (i in 1:repl.col) {
			if (debuglevel>=5) 
				cat("\nColumn:", i, ", method: transect\n--------------------")
			res <- do.call(paste('fd.estim.', method1d, sep=''), c(list(data[,i], p.index=p), args))
			fd[,i+repl.row] <- c(res$fd, res$scale)
		}
	}
	if(plot.loglog)
		warning('No plotting feature for the method transect')

	DD <- list(	fd = 1 + median(fd['fd',]),
	  			scale=median(fd['scale',]))
	  			
	invisible(createFractalDim(dim=1, methods=method.name, FD=DD,
                              window.size=repl.row, data.dim=2))
}

fd.estim.transect.var <- function(data, p.index=2, ...) {
	return (do.fd.estimate2d.simple(data, p=p.index, method.name='transect.var', method1d='variation', ...))
}

fd.estim.transect.incr1 <- function(data, p.index=2, ...) {
	return (do.fd.estimate2d.simple(data, p=p.index, method.name='transect.incr1', method1d='incr1',...))
}

fd.estim.isotropic <- function(data, p.index=2, direction='hvd+d-',
							plot.loglog=FALSE, nlags = "auto", 
							plot.allpoints=FALSE, legend.type='s', ...,
                          	debuglevel=0) {
    # Compute the isotropic empirical variogram (or the equivalent for the variation estimator),  
	# by averaging over given directions (horizontal, vertical, diagonal ...). 
	# Then apply log-log regression. See page 12 in Davies and Hall (1999).
	return(do.fd.estimate2d.method(data, 'fd.Variation.isotropic', method.name='isotropic',
										method.name.main='Isotropic', p=p.index, direction=direction,
										plot.loglog=plot.loglog, nlags=nlags,
										plot.allpoints=plot.allpoints,
										legend.type=legend.type, ..., debuglevel=debuglevel))

}

fd.estim.squareincr <- function(data, p.index=2,
							plot.loglog=FALSE, nlags = "auto", 
							plot.allpoints=FALSE, legend.type='s', ...,
                          	debuglevel=0) {
	return(do.fd.estimate2d.method(data, 'fd.Variation.squareincr', method.name='squareincr',
										method.name.main='Square Increment', 
										p=p.index, direction='h',
										plot.loglog=plot.loglog, nlags=nlags,
										plot.allpoints=plot.allpoints,
										legend.type=legend.type, ..., debuglevel=debuglevel))

}

fd.estim.filter1 <- function(data, p.index=2, direction='hvd+d-',
							plot.loglog=FALSE, nlags = "auto", 
							plot.allpoints=FALSE, legend.type='s', ...,
                          	debuglevel=0) {
	return(do.fd.estimate2d.method(data, 'fd.Variation.filter1', method.name='filter1', 
										method.name.main='Filter1', 
										p=p.index, direction=direction,
										plot.loglog=plot.loglog, nlags=nlags,
										plot.allpoints=plot.allpoints, 
										legend.type=legend.type, ..., debuglevel=debuglevel))

}


do.fd.estimate2d.method <- function(data, fd.function, method.name, p=2, direction='hvd+d-',
							method.name.main = '',
							plot.loglog=FALSE, nlags = "auto", 
							plot.allpoints=FALSE, legend.type='s', ...,
                          	debuglevel=0) {
                          		
	n<-attr(data,"dim")[1]
 	if (n != attr(data,"dim")[2]) 
    	stop("The method works with square matrices only.")
  	x<-1:(n-1)
  	diagonal <- has.direction(direction, 'd')
  	horver <- has.direction(direction, 'v') || has.direction(direction, 'h')
  	if (diagonal) { # is 'diagonal' one of the directions,
    		lagsm<-x*sqrt(2)         # put multiple of sqrt(2) on the x-axis
    		if (!horver) x<-lagsm
    		else x<-sort(c(x,lagsm))
  	}
  	lx <- length(x)
	all.x <- x
	
  	if (nlags == "auto") {
  		if (diagonal && horver) nr.lags <- 3
  		else nr.lags <- 2	
  	} else {
  		if (nlags == "all")
  			nr.lags <-  lx-1 
  		else {
  			nlags <- as.numeric(nlags)
  			nr.lags <- min(length(x[x<= nlags]), lx-1)
  		}
  	}
  	
	x<-x[1:nr.lags] 	
    g <- do.call(fd.function, list(data=data, lags=x, p=p, direction=direction))
     	
    logx<-log(x)
  	logg <- log(g)
  	
    rawD <- get.rawFD.from.regression(logx, logg)
  	
	if (debuglevel>=5) {
    	cat("\nLog-log plot output:\n")
    	summary(rawD)
  	}
  	if (plot.loglog) {
  		grest <- NULL
    	if (plot.allpoints & lx-1 > nr.lags) {
      		grest <- rep(NA, lx-1-nr.lags)
      		for (ll in 1:(lx-1-nr.lags)) {
        		grest[ll] <- do.call(fd.function, list(data=data, lags=all.x[nr.lags+ll], 
        								p=p, direction=direction))
      		}
      		logx <- log(all.x[1:(lx-1)])
      		logg <- c(log(g), log(grest))
      		isnotna <- !is.nan(logg)
      		rawD$x <- logx[isnotna]
      		rawD$y <- logg[isnotna]
    	}
    	plot.args <- list(...)
    	if (!is.element('main', names(plot.args))) plot.args[['main']] <- paste(method.name.main, "- p =",p)
    	if (!is.element('ylab', names(plot.args))) plot.args[['ylab']] <- "log(p-variation)"
    	if (!is.element('xlab', names(plot.args))) plot.args[['xlab']] <- paste("log(lags = [1:",
                         round(exp(max(rawD$x)),2),"])", sep="")
        filled <- if(plot.allpoints && !is.null(grest)) 1:nr.lags else NULL
    	do.call('fd.plot.loglog', c(list(rawD=rawD, filled=filled), plot.args))
  	}
  	DD <- list(	fd = if (!is.na(rawD$alpha)) 3-rawD$alpha/p else NA,
  				# normalized intercept
  				scale = if(!is.na(rawD$intercept)) exp(rawD$intercept)^(1/p) else NA				)
  
  	if (plot.loglog) fd.plot.legend(rawD, DD, legend.type=legend.type)
  	
  	invisible (createFractalDim(dim=1, methods=method.name, FD=DD,
                              window.size=n, data.dim=2, loglog=rawD))

}


fd.Variation.isotropic <- function(data, p=2, lags=1:2,  direction='h') {
	lx <- attr(data,"dim")[1]
  	ly <- attr(data,"dim")[2]
  	m <- length(lags)
  	g<-rep(0, m)
  	int.lags <- 1:max(lags)
  	for (j in 1:m) {
  		count <- 0
  		if(is.element(lags[j], int.lags)) { # lag is an integer
  			if (has.direction(direction, 'v')) {  # vertical
      			g[j]<- g[j] + sum(abs(data[(lags[j]+1):lx,1:ly]-data[1:(lx-lags[j]),1:ly])^p)/((lx-lags[j])*ly)
      				count <- count+1
   			}
   			if (has.direction(direction, 'h')) {  # horizontal
  				g[j]<- g[j] + sum(abs(data[1:lx,(lags[j]+1):ly]-data[1:lx,1:(ly-lags[j])])^p)/((ly-lags[j])*lx)
  				count <- count+1
   			}
   		} else {
   			for(dir in c('d+', 'd-')) {
   				if (has.direction(direction, dir)) {
   					g[j]<- g[j] + fd.VariationFiltered(data, p=p, lags=lags[j], direction=dir)
   					count <- count+1
   				}
   			}
   		}
   		if (count > 0)
			g[j] <- g[j]/count
	}
  return(g)
}

fd.Variation.squareincr <- function(data, p=2, lags=1:2, ...) {
	return(fd.VariationFiltered(data, p=p, lags=lags, direction='h', filter=3))
}

fd.Variation.filter1 <- function(data, p=2, lags=1:2,  direction='h') {
	return(fd.VariationFiltered(data, p=p, lags=lags, direction=direction, filter=1))
}

fd.VariationFiltered <- function(data, p=2, lags=1:2, direction='hvd+d-', filter = 0) {
	# Filtering for increments of order 0, according to 
  	# Zhu, Stein: Parameter Estimation for Fractional
  	# Brownian Surfaces, Statistica Sinica 12, 2002.
  
  	lx <- attr(data,"dim")[1]
  	ly <- attr(data,"dim")[2]
  	m <- length(lags)
  	seq.int<-1:max(lags)
  	trans <- matrix(0,nrow=4, ncol=4)
  	trans[1,]<-c(1,0,0,1)  # I
  	trans[2,]<-c(0,1,-1,0) # R 90' rotation
  	trans[3,]<-c(1,1,-1,1) # D 45' rotation + scaling
  	trans[4,]<-c(-1,1,1,1) # 135' rotation
  	Tdirs <- c('h', 'v', 'd+', 'd-')
  	filters <- list()
  	filters[['0']] <- list(A=c(1,-1), delta=matrix(c(1,0,0,0), 2,2), M=2)
  	filters[['1']] <- list(A=c(1,1,-2), delta=matrix(c(1,-1,0,0,0,0), 3,2), M=2)
  	filters[['3']] <- list(A=c(1,1,-1,-1), delta=matrix(c(1,0,1,0,0,1,1,0), 4,2), M=1)
  
  	use.filter <- filters[[as.character(filter)]]
	delta <- use.filter$delta
	ndelta <- nrow(delta)
	A <- use.filter$A
	lf <-length(A)
	M <- use.filter$M
	
  	g<-rep(0,m)
  	for (ilag in 1:m) {
    	lag<-lags[ilag]
   	 	tmp<-0
    	if (is.element(lag, seq.int)) { # lag is an integer
      		denom <- 1
      		tidx <- 0
    	} else {  # lag is real (diagonal rotation)
      		tidx <- M
      		denom <- sqrt(2)
    	}
    	iM<-0
    	for (tran in 1:M) {
      		if (!has.direction(direction, Tdirs[tran+tidx])) next
        	iM<-iM+1
        	deltatran <- matrix(NA, nrow=2, ncol=ndelta)
        	for (j in 1:ndelta) {
        		deltatran[,j] <-  t((lag/denom*matrix(trans[tran+tidx,],2,2)) %*% delta[j,])
        	}
			# shift transformed points to the area
            # of positive coordinates
            for(i in 1:nrow(deltatran)) {
          		mind <- min(deltatran[i,])
          		if (mind < 0) {
            		deltatran[i,] <- deltatran[i,] - mind
          		}
          	}
          	this.result <- 0
        	maxx<-max(deltatran[1,])
        	maxy<-max(deltatran[2,])
			dt2 <- 1 + deltatran
			xd<-lx-maxx+deltatran
			yd<-ly-maxy+deltatran
			if (max(dt2) > lx) {
          		g[ilag]<-NaN
          		break
        	}
        	for (j in 1:ndelta) {
        		this.result <- this.result + A[j]*data[dt2[1,j]:xd[1,j],dt2[2,j]:yd[2,j]]
            }
            tmp <- tmp + sum(abs(this.result)^p)/((lx-maxx)*(ly-maxy))
            
    	}
    	if (iM==0) g[ilag]<-NaN
    	else {
      		if (!is.nan(g[ilag]))
        		g[ilag]<-tmp/iM
    	}
  	}
  	return(g)   
}







