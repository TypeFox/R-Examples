fd.get.available.methods <- function(dim=1) {
	switch(dim, 
		list( 'variogram', 'madogram', 'rodogram', 'variation', 'incr1', 'boxcount', 'hallwood', 'periodogram', 'wavelet',
             'dctII', 'genton'),
        list('transect.var', 'transect.incr1', 'isotropic', 'squareincr', 'filter1')
        )
}

fd.get.method.coding <- function(dim=1) {
	methods <- fd.get.available.methods(dim=dim)
  	return (1:length(methods))
}
  
get.vector.method.coding <- function(methods, dim=1) {
	# get codes for given methods. 
  	coding <- fd.get.method.coding(dim=dim)
  	return(coding[is.element(fd.get.available.methods(dim), methods)])
}

get.methods.from.coding <- function(coding, dim=1) {
	return(unlist(fd.get.available.methods(dim)[coding]))
}

"fd.estimate" <- function (data, ...) UseMethod("fd.estimate")

fd.estimate.numeric <- function (data, methods="madogram", 
							window.size=length(data),
                   			step.size=window.size, trim=TRUE, 
                   			keep.data=FALSE, keep.loglog=FALSE,
                   			parallel = FALSE, nr.nodes=NULL, debuglevel=0, ...) {  
  	# 1d data
    args <- prepare.for.row.estimation(data, methods, window.size,
                                       		step.size, dim=1, ...)
    if (debuglevel>=1)
      	cat("\nNumber of iterations: ", args$niterations,"\n")
    	
    fractaldim <- do.fd.estimate(args=args, data=data, methods=methods,
                           step.size=step.size, trim=trim, keep.loglog=keep.loglog, 
                           parallel=parallel, nr.nodes=nr.nodes,
                           debuglevel=debuglevel, ...)
                           
    dim <- if (args$niterations==1) 1 else 2
    datatmp <- if (keep.data) data else NULL
    return(invisible(structure(list(dim=dim, fd=fractaldim$fd, 
                             scale=fractaldim$scale,
                             methods = methods,
                             methods.coding=args$method.coding, data=datatmp,
                             data.dim = 1,
                             window.size=args$window.size,
                             step.size=step.size,
                             loglog=fractaldim$loglog),
                        class="FractalDim")))
}

"fd.estimate.data.frame" <- function (data, ...) UseMethod("fd.estimate.matrix")

fd.estimate.matrix <- function (data, methods="transect.var",
							window.size=ncol(data),
                   			step.size=window.size, trim=TRUE, keep.data=FALSE,
                   			keep.loglog=FALSE, parallel = FALSE, nr.nodes=NULL, 
                   			debuglevel=0, ...) {
	# 2d data
    repl <- trunc((nrow(data)-window.size)/step.size)+1
  	rephor <- trunc((ncol(data)-window.size)/step.size)+1 # (horizontal along the columns)

    if (debuglevel>=1) {
    	cat(paste("Total iterations over rows:" , repl,"\n"))
    	cat(paste("Total iterations over columns:" , rephor,"\n"))
    }
    
    args <- prepare.for.row.estimation(data[1,], methods, window.size,
                                       step.size, dim=2, ...)
    complete.result <- list()
    complete.result[['loglog']] <- list()
    result.names <- c('fd', 'scale')
    if (!parallel) { # run sequentially (one row after another)
    	for (it in 1:repl) {
        	if (debuglevel>=2)
          		cat(paste("Iteration:" , it,"\n"))
        	fractaldim <- do.fd.estimate2d( it, args=args, data=data,
                                     step.size=step.size, trim=trim,
                                     keep.loglog=keep.loglog, debuglevel=debuglevel)
          	for (what in result.names) {
           		complete.result[[what]] <- abind(complete.result[[what]],
                                           		array(fractaldim[[what]], 
                                           		c(1,dim(fractaldim[[what]]))), 
                                           	along=1)
          	}
          	if(keep.loglog)
          		complete.result$loglog[[it]] <- fractaldim$loglog
        }
    } else { # run rows in parallel
    	require(snowFT)
    	nr.nodes <- if(is.null(nr.nodes)) repl else min(repl, nr.nodes)
      	fractaldims <- performParallel(nr.nodes, 1:repl,
                                     do.fd.estimate2d,
                                     args=args, data=data,  methods=methods,
                                     step.size=step.size, trim=trim, 
                                     keep.loglog=keep.loglog,
                                     debuglevel=debuglevel, gentype='None',
                                     initfun = fd.initfunction)
      	for (what in result.names) {
        	complete.result[[what]] <- array(fractaldims[[1]][[what]],
                                         c(1, dim(fractaldims[[1]][[what]])))
      	}
      	if(keep.loglog)
          	complete.result$loglog[[1]] <- fractaldim[[1]]$loglog
      	for (it in 2:repl) {
        	for (what in result.names) {
          		complete.result[[what]] <- abind(complete.result[[what]],
                                           fractaldims[[it]][[what]], along=1)
        	}
        	if(keep.loglog)
          		complete.result$loglog[[it]] <- fractaldim[[it]]$loglog
      	}
    }
    if (debuglevel >= 2)
    	cat("Iterating over rows done.\n")
    datatmp <- if (keep.data) data else NULL

    invisible(structure(list(dim=3, fd=complete.result$fd,
                             scale=complete.result$scale,
                             methods = methods,
                             methods.coding=args$method.coding, data=datatmp,
                             data.dim = 2,
                             window.size=args$window.size,
                             step.size=step.size,
                             loglog = complete.result$loglog),
                        class="FractalDim"))
}

				
prepare.for.row.estimation <- function (data, methods, window.size, step.size, dim=1, ...) {
  	n <- length(data)
  	nest<-length(methods)
  	if (nest <= 0) {
    	stop("No estimation method specified.")
  	}
  	window.size <- max(2, min(n, window.size))
  	niterations <- trunc((n - window.size)/step.size) + 1
  	optargs <- list(...)
  	methodlist <- list()
  	method.coding <- c()
  	for (imeth in 1:nest) {
  		this.meth <- methods[[imeth]]
  		args <- NULL
  		if (is.list(this.meth)) {
  			name <- this.meth[[1]]
  			if(length(this.meth) > 1) args <- this.meth[2:length(this.meth)]
  		} else name <- this.meth
      	methodlist <- c(methodlist, list(list(name=name, args=c(args, optargs))))
      	method.coding <- c(method.coding, get.vector.method.coding(name, dim=dim))
  	}
  	return(list(n=n, nest=nest, window.size=window.size, niterations=niterations,
              methodlist=methodlist, method.coding=method.coding))
}

do.fd.estimate.for.parallel <- function (repl, args, data,
                                         methods=c("variation"),
                                         step.size=args$window.size, trim=TRUE,
                                         debuglevel=0, ...) {
  return(do.fd.estimate( args=args, data=data[repl,],
                        methods=methods,
                        step.size=step.size, trim=trim, 
                        parallel=FALSE, # since this function is already running in parallel
                        debuglevel=debuglevel, ...))
}

do.fd.estimate <- function (args, data, methods=c("variation"),
                            step.size=args$window.size, trim=TRUE, keep.loglog=FALSE,
                            parallel=FALSE, nr.nodes=NULL, debuglevel=0, ...) {

	niterations <- args$niterations
  	nest <- args$nest
  	methodlist <- args$methodlist
  	window.size <- args$window.size
  	if(parallel) {
  		require(snowFT)
  		nr.nodes <- if(is.null(nr.nodes)) niterations else min(niterations, nr.nodes)
  		resultlist <- performParallel(nr.nodes, 1:niterations, do.fd.estimate.1iteration, 
  									data=data, step.size=step.size, 
  									window.size=window.size, nest=nest, 
  									methodlist=methodlist, trim=trim, keep.loglog=keep.loglog, 
  									debuglevel=debuglevel, gentype='None',
                                    initfun = fd.initfunction, ...)	
  	} else { # sequentiell
  		resultlist <- list()
  		for (sl in 1:niterations) {  # moving sliding window 
  			resultlist[[sl]] <- do.fd.estimate.1iteration(sl, data=data, step.size= step.size, 
  										window.size=window.size, nest=nest, 
  										methodlist=methodlist, trim=trim, 
  										keep.loglog=keep.loglog, debuglevel=debuglevel)
  		}
  	}
  	result <- list(fd = matrix(0,nrow=niterations,ncol=nest),
  				scale = matrix(0,nrow=niterations,ncol=nest),
  				loglog = list(),
  				method.args = methodlist
  			)
  	for (sl in 1:niterations) {
		result$fd[sl,]<-resultlist[[sl]]$fd
        result$scale[sl,]<-resultlist[[sl]]$scale
        if(keep.loglog)
        	result$loglog[[sl]] <- resultlist[[sl]]$loglog
     }
  return (result)
}

do.fd.estimate.1iteration <- function(it, data, step.size, window.size, nest,
										methodlist, trim, keep.loglog, debuglevel=0) {
	if (debuglevel>=3)
    	cat("window nr. ", it,"\n")
	result <- list(fd=rep(0, nest), scale=rep(0, nest), loglog=list())
    dat <- data[((it-1)*step.size+1):((it-1)*step.size+window.size)]
    if (length(dat[is.na(dat)]) > 0) return(result) # don't estimate if NAs present
    for (imeth in 1:nest) {
    	this.meth <- methodlist[[imeth]]$name
        if (debuglevel >= 4)
          	cat("Method: ", this.meth,"\n")
        est <- do.call(paste("fd.estim.", this.meth, sep=""),
                         c(list(data=dat, debuglevel=debuglevel), methodlist[[imeth]]$args))
        result$fd[imeth]<-est$fd
        result$scale[imeth]<-est$scale
        if(keep.loglog)
        	result$loglog[[imeth]] <- est$loglog
	}
    if(trim) {
		# trim (estimate in [1,2])
    	result$fd<-pmin(pmax(1,result$fd),2) 
	}
	return(result)
}
	
do.fd.estimate2d <- function (iteration, args, data, 
                            step.size=args$window.size, trim=TRUE, 
                            keep.loglog=FALSE, debuglevel=0) {

	niterations <- args$niterations
  	nest <- args$nest
  	methodlist <- args$methodlist
  	window.size <- args$window.size
  	window.size.vertical <- min(window.size, dim(data)[1]) 
    data <- data[((iteration-1)*step.size+1):((iteration-1)*step.size+window.size.vertical),]
  	columns <- ncol(data)

  	resultlist <- list()
  	for (sl in 1:niterations) {  # moving sliding window horizontally
  		if (debuglevel>=3) 
      		cat("\nwindow nr. ", sl)
		result <- list(fd=rep(0, nest), scale=rep(0, nest), loglog=list())
    	dat <- data[,((sl-1)*step.size+1):((sl-1)*step.size+window.size)]
    	if (length(dat[is.na(dat)])>0) next
    	for (imeth in 1:nest) {
    		this.meth <- methodlist[[imeth]]$name
        	if (debuglevel >= 4)
          		cat("\nMethod: ", this.meth)
          	est <- do.call(paste("fd.estim.", this.meth, sep=""),
                         c(list(data=dat, debuglevel=debuglevel), methodlist[[imeth]]$args))
        	
        	result$fd[imeth]<-est$fd
        	result$scale[imeth]<-est$scale
        	if(keep.loglog)
        		result$loglog[[imeth]] <- est$loglog
		}
		if (debuglevel >= 4)
          	cat("\n")
      	if(trim) {
			# trim (estimate in [2,3])
      		result$fd<-pmin(pmax(2,result$fd),3) 
      	}
		resultlist[[sl]] <- result
    }
  	result <- list(fd = matrix(0,nrow=niterations,ncol=nest),
  				scale = matrix(0,nrow=niterations,ncol=nest),
  				loglog = list(),
  				method.args = methodlist
  			)
  	for (sl in 1:niterations) {
		result$fd[sl,]<-resultlist[[sl]]$fd
        result$scale[sl,]<-resultlist[[sl]]$scale
        if(keep.loglog)
        	result$loglog[[sl]] <- resultlist[[sl]]$loglog
     }
  return (result)
}
	
	
createFractalDim <- function (dim, methods, FD,
                              methods.coding=get.vector.method.coding(methods),
                              window.size=0,
                              step.size=window.size, data.dim=1, loglog=NULL) {
  	return (structure(list(dim=dim, 
  						fd=FD$fd, 
  						scale=FD$scale, 
                        methods=methods, methods.coding=methods.coding,
                        window.size=window.size,
                        step.size=step.size,
                        data.dim=data.dim,
                        loglog=loglog), class="FractalDim"))


}


create.plot.legend <- function(D, rawD=NULL, legend.type='s') {
	if(legend.type=='n') return(NULL)
	if(legend.type=='f') # full legend
    	label <- substitute(paste(hat(D) == d,  '(', al, '), ', scale==s, '(', inter, ')', sep=''),
                   list(d=round(D$fd,2), al=round(rawD$alpha,2),
                   		s=round(D$scale,2), inter=round(rawD$intercept,2)))	else {
    	if(legend.type=='s')          			
    		label <- substitute(hat(D) == d, list(d=format(round(D$fd,2), nsmall=2)))
    	else stop('Unknown value of legend.type. Allowed: "f" (full), "s" (short), "n" (None)')
    }
	return(label)
} 


fd.regression <- function (x, y, leaveout=0) {
  xx <- x[(leaveout+1):length(x)]
  yy <- y[(leaveout+1):length(y)]
  xbar <- mean(xx)
  slope <- sum((xx-xbar)*yy)/sum((xx-xbar)^2)
  ybar <- mean(yy)
  intercept <- ybar-slope*xbar
  lsq <- sum((yy-(intercept+slope*xx))^2)/length(xx)
  return (list(alpha=slope, intercept=intercept, lsq=lsq))
}

get.rawFD.from.regression <- function (x, y, leaveout=0) {
  # Returns an object FDloglog with values obtained from
  # a linear regression
  regression <- fd.regression(x, y, leaveout)
  return (structure(list(alpha=regression$alpha,
                         intercept=regression$intercept,
                         x=x, y=y, n=length(x),
                         lsq=regression$lsq),
                    class="FDloglog"))
}
      
  
fd.initfunction <- function() {
  library(fractaldim)
}

fd.get <- function(fractaldim, method) {
	beta <- get.vector.method.coding(c(method), dim=fractaldim$data.dim)
	method.index <- (1:length(fractaldim$methods.coding))[fractaldim$methods.coding==beta]
	if (length(method.index)<=0) {
		cat("\nEstimate for method ",method," not found.\n")
		return(NULL)
	}
	get.index2d <- function(x) {
		return(sapply(method.index, function(y) c(x, y)))
	}
	get.index3d <- function(x) {
		do.get.index3d <- function(y) return(sapply(method.index, function(z) c(x, y, z)))
		return(sapply(1:dim(fractaldim$fd)[2], do.get.index3d))
	}
	index <- switch(fractaldim$dim,  
  					method.index,   # dim=1 
  					t(sapply(1:nrow(fractaldim$fd), get.index2d)),  # dim=2
  					t(matrix(sapply(1:dim(fractaldim$fd)[1], get.index3d), nrow=3)))
  	loglog <- NULL
	if (!is.null(fractaldim$loglog)) loglog <- lapply(fractaldim$loglog, function(x) x[method.index])
  	
	result <- structure(list(dim=fractaldim$dim, methods=method,
  								methods.coding=beta,
  								window.size=fractaldim$window.size,
  								step.size=fractaldim$step.size, 
  								fd=fractaldim$fd[index],
  								scale=fractaldim$scale[index],
  								data.dim=fractaldim$data.dim,
  								data=fractaldim$data,
  								loglog=loglog,
  								method.args=fractaldim$method.args[method.index]
  								), class='FractalDim')
  								
	if (fractaldim$dim >= 2) {
		d <- dim(fractaldim$fd)
		ld <- length(dim(fractaldim$fd))
		# replace the last dimension by the number of new methods
		for (what in c("fd","scale")) {
	   		result[[what]] <- array(result[[what]], c(d[-ld], length(method.index)))
		}
	}
  	return(result)							
}



fd.plot.loglog <- function(rawD, col="red", filled=NULL, ...) {
	x <- rawD$x
	y <- rawD$y
	plot(x=x, y = y, ...)
	if (!is.null(filled)) 
		points(x[filled], y[filled], bg='black', pch=21)
	abline(rawD$intercept, rawD$alpha, col=col)
}

fd.plot.legend <- function(rawD, D, legend.type='s') {
	legend.text <- create.plot.legend(D, rawD, legend.type=legend.type)
  	if(!is.null(legend.text)) {
  		where <- if(rawD$alpha > 0) 'topleft' else 'topright' # gradient direction
  		legend(where, legend=legend.text, bty='n')
  	}
}

summary.FDloglog <- function (object, ...){
  output <- matrix(c(object$alpha, object$intercept, object$lsq,
                     object$n), 1,4,
      dimnames=list(c("regression"),
        c("slope", "intercept", "least square", "size")))
  print(output)
  invisible()
}

summary.FractalDim <- function (object, ...){
  s <- "Grid: "
  if (object$dim == 1) {
    s<-paste(s,1)
    means.fd <- object$fd
    means.int <- object$scale
    sd.fd <- sd.int <- rep(NA, length(object$fd))
  } else {
    s<-paste(s,dim(object$fd)[1])
    if (object$dim > 2) {
      s<-paste(s, "x", dim(object$fd)[2])
    }
    along <- object$dim
    means.fd <- apply(object$fd, c(along), mean, na.rm=TRUE)
    means.int <- apply(object$scale, c(along), mean, na.rm=TRUE)
    sd.fd <-  sd.int <- rep(NA,length(means.fd))
    if (dim(object$fd)[1] > 1) {
      if (object$dim == 3) { # flatten 1. and 2. dim,
                           # otherwise sd returns a matrix      
        tmp <- apply(object$fd, c(along), c)
        sd.fd <- apply(tmp, c(2), sd, na.rm=TRUE)
        tmp <- apply(object$scale, c(along), c)
        sd.int <- apply(tmp, c(2), sd, na.rm=TRUE)
      } else {
        sd.fd <- apply(object$fd, c(along), sd, na.rm=TRUE)
        sd.int <- apply(object$scale, c(along), sd, na.rm=TRUE)
      }
    }
  }
  cat("\n",s,"\n\n")
	method.names <- c()
	for (meth in object$methods) {
		if (is.list(meth)) {
			this.meth <- meth$name
			for (arg in names(meth)) {
				if (arg != 'name') this.meth <- paste(this.meth, ', ', 
									arg, '=', paste(meth[[arg]], collapse=','), sep='')
			}
			method.names <- c(method.names, this.meth)
		} else method.names <- c(method.names, meth)
	}
  r<-3
  output <- matrix(c(round(means.fd,r),
                      round(sd.fd,r), 
                      round(means.int,r),
                      round(sd.int,r)),
                    length(object$methods.coding), 4,
                    dimnames=list(method.names, c("mean(fd)", "sd(fd)",
                      "mean(sc)", "sd(sc)")))
  print(output, quote=FALSE)
  invisible()
}

dim.FractalDim <- function(x) {
  return (x$dim)
}
