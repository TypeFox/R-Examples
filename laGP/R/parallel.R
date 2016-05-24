#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************


## agp.chunk:
##
## a version of the aGP function which works with a subset of the XX
## locations, for use with aGP.parallel for parallelization

aGP.chunk <- function(chunk, XX, X, Z, start, end, d, g, method, 
	Xi.ret, close, numrays, num.gpus, gpu.threads, omp.threads, nn.gpu, verb)
  {
    ## might need to chunk-up the d argument
    if(!is.null(d)) {
      if(!is.list(d) && length(d) > 1) d <- d[chunk]
      else if(is.list(d) && !is.null(length(d$start)) && length(d$start) > 1) 
        d$start <- d$start[chunk]
    }

    ## might need to chunk up the g argument
    if(!is.null(g)) {
      if(!is.list(g) && length(g) > 1) g <- g[chunk]
      else if(is.list(g) && !is.null(g$start) && length(g$start > 1)) 
        g$start <- g$start[chunk]
    }

    ## might need to chunk up nn.gpu
    if(nn.gpu == nrow(XX)) nn.gpu <- length(chunk)
    else {
      p <- nn.gpu / nrow(XX)
      nn.gpu <- round(length(chunk) * p)
    }

    ## run with chunked up XX
    aGP(X, Z, XX=XX[chunk,], start, end, d, g, method, Xi.ret, 
    	close, numrays, num.gpus, gpu.threads, omp.threads, nn.gpu, verb)
  }


## aGP.parallel:
##
## a version of the aGP function for use with clusters made by snow
## and/or the parallel package.  Uses clusterApply to split up the
## XX rows into a particular number of chunks, and then combines
## the results into a single aGP output object

aGP.parallel <- function(cls, XX, chunks=length(cls), X, Z, start=6, end=50, d=NULL, 
                     g=1/1000, method=c("alc", "alcray", "mspe", "nn", "efi"), 
                     Xi.ret=TRUE, close=min(1000*if(method == "alcray") 10 else 1, nrow(X)),
                     numrays=ncol(X), num.gpus=0, gpu.threads=num.gpus, 
                     omp.threads=if(num.gpus > 0) 0 else 1, 
                     nn.gpu=if(num.gpus > 0) nrow(XX) else 0, verb=1)
  {
    ## timing
    tic <- proc.time()[3]

    ## creates length(cls) chunks
    ## all.chunks <- clusterSplit(cls, 1:n) 
    ## creates user specified number of chunks
    n <- nrow(XX)
    all.chunks <- split(1:n, rep(1:chunks, length=n))

    ## compute in parallel
    clusterEvalQ(cls, library(laGP))
    all.outs <- clusterApply(cls, all.chunks, aGP.chunk, XX=XX, X=X, Z=Z, 
    	start=start, end=end, d=d, g=g, method=method, Xi.ret=Xi.ret, 
    	close=close, numrays=numrays, num.gpus=num.gpus, gpu.threads=gpu.threads,
     	omp.threads=omp.threads, nn.gpu=nn.gpu, verb=verb)

    ## allocate a single output object
    out <- list(
    	mean=rep(NA, n),
    	var=rep(NA, n),
    	llik=rep(NA, n),
    	d=all.outs[[1]]$d,
    	g=all.outs[[1]]$g,
    	chunk.times=rep(NA,chunks),
    	method=all.outs[[1]]$method)

    ## allocate optional outputs for single output object
    ## mle
    if(!is.null(all.outs[[1]]$mle)) {
    	out$mle <- matrix(NA, nrow=nrow(XX), ncol=ncol(all.outs[[1]]$mle))
    	out$mle <- as.data.frame(out$mle)
    	names(out$mle) <- names(all.outs[[1]]$mle)
    }
    ## Xi indices
    if(!is.null(all.outs[[1]]$Xi))
    	out$Xi <- matrix(NA, nrow=nrow(XX), ncol=ncol(all.outs[[1]]$Xi))

    ## loop over each output and added data into the combined output object
    for(i in 1:length(all.outs)) {
    	out$mean[all.chunks[[i]]] <- all.outs[[i]]$mean
     	out$var[all.chunks[[i]]] <- all.outs[[i]]$var
     	if(!is.null(out$mle)) 
     		out$mle[all.chunks[[i]],] <- all.outs[[i]]$mle
     	if(!is.null(out$Xi)) out$Xi[all.chunks[[i]],] <- all.outs[[i]]$Xi
     	out$chunk.times[i] <- all.outs[[i]]$time
    }

    ## done timing
    toc <- proc.time()[3]
    out$time <- toc - tic

    ## all done
	return(out)
  }
