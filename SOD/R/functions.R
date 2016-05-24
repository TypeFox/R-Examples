##
## A set of functions that are useful for visualising squeezed data ##

##loadModule("mod_R_DimSqueezer", TRUE)
## Module is dangerous; don't load, but
## use the following functions instead

## Use reference classes:
DimSqueezer <- setRefClass("DimSqueezer",
                           fields=list(
                               pointer="externalptr",
                               data.matrix="matrix"))
DimSqueezer$methods(initialize =
                    function(m){
                        data.matrix <<- m
                        pointer <<- .Call("DimSqueezer", m, "SOD")
                    },
                    squeeze =
                    function(target_dim, iter_no){
                        .Call("squeeze", pointer, target_dim, iter_no, "SOD")
                    },
                    squeezeDF =
                    function(dimFactors){
                        .Call("squeezeDF", pointer, dimFactors, "SOD")
                    },
                    useOpenMP =
                    function(useOMP){
                        invisible(.Call("useOpenMP", pointer, useOMP, "SOD"))
                        
                    },
                    setThreadNo =
                    function(threadNo){
                        invisible(.Call("setThreadNo", pointer, threadNo, "SOD"))
                    },
                    residual =
                    function(remResidual){
                        .Call("removeResidualStress", pointer, remResidual, "SOD")
                    }
                    )
DimSqueezer$lock("pointer")

## note that we don't need to register a finalizer as the
## Rcpp ExtPtr class seems to take care of that.
if(TRUE){
    ## Repeat for the CL based Class
    DimSqueezer_CL <- setRefClass("DimSqueezer_CL",
                                  fields=list(
                                      pointer="externalptr",
                                      data.matrix="matrix"))
    
    DimSqueezer_CL$methods(initialize =
                           function(m){
                               data.matrix <<- m
                               pointer <<- .Call("DimSqueezer_CL", m, "SOD")
                           },
                           squeeze =
                           function(target_dim, iter_no, wksize){
                               .Call("squeeze_cl", pointer, target_dim, iter_no, wksize)
                           }
                           )
    
    ## note that we don't need to register a finalizer as the
    ## Rcpp ExtPtr class seems to take care of that.
    DimSqueezer_CL$lock("pointer")
}

## this is good for safety, but does not allow me to automatically
## rebuild the pointer from the arguments.

plot.sod_sq3 <- function(x, ..., ptype='p'){
    p.type = NULL
    mc <- match.call()
    switch(ptype,
           s = {p.type = "plotStress"},
           c = {p.type = "plotConcentric"},
           p = {p.type = "plotPoints"},
           warning("Unknown plot option")
       )
    if(!is.null(p.type)){
        mc[[1]] <- as.name(p.type)
        mc$ptype = NULL
        eval(mc)
    }
}

## makes a color for each of level of v
## with low (blue) to high (purple) via, cyan, green, yellow, red.
## this can also be done by reordering the
## output of the rainbow function.
## but not sure how to get the radial shift.
hsvScale <- function(v, sat=1, val=0.75, alpha=1, min.v=min(v), max.v=max(v)){
  ## run from blue (4/6) -> magenta (5/6)
  v.range <- max.v - min.v
  if(!v.range)
    return(rep(hsv(4/6), length(v)))
  v <- 5 + 5 * (min.v - v) / v.range
  ## now runs from 5 (magenta) -> 0 (red)
  ## convert to 4, 3, 2, 1, 0, 5
  v <- (v - 1) ## and now runs 4, 3, 2, 1, 0, -1
  v[ v < 0 ] <- 6 + v[ v < 0 ] ## -> 4, 3, 2, 1, 0, 5
  hsv( v/6, sat, val, alpha )
}


## sq squeezed data
plotPoints <- function(x, col=hsvScale(x$node_stress), xc=1, yc=2,
                        invert.y=FALSE, xlab=NA, ylab=NA, ...){
    xv = x$pos[,xc]
    yv = x$pos[,yc]
    if(invert.y)
        yv = -yv
    plot(xv, yv, col=col, bg=col, xlab=xlab, ylab=ylab, ...)
}

plotConcentric <- function(x, cex.data, col=hsvScale(1:ncol(cex.data)),
                            xc=1, yc=2, cex.max=3, invert.y=FALSE, pch=19, xlab=NA, ylab=NA, leg.pos=NULL, ...){
    xv <- x$pos[,xc]
    yv <- x$pos[,yc]
    if(invert.y)
        yv <- -yv
    p.cex <- matrix(nrow=nrow(cex.data), ncol=ncol(cex.data), data=0)
    p.cex[,ncol(cex.data)] <- sqrt(cex.data[,ncol(cex.data)])
    for(i in (ncol(cex.data)-1):1){
        p.cex[,i] <- p.cex[,(i+1)] + sqrt(cex.data[,i])
    }
    ## scale to cex.max
    p.cex <- cex.max * p.cex / max(p.cex)
    plot(xv, yv, type='n', xlab=xlab, ylab=ylab, ...)
    for(i in 1:ncol(p.cex))
        points(xv, yv, cex=p.cex[,i], col=col[i], bg=col[i], pch=pch)

    if(!is.null(leg.pos)){
        legend(leg.pos, legend=colnames(cex.data), col=col, pch=pch)
    }
}

plotStress <- function(x, bg.alpha=0.5, bg.sat=1, bg.val=0.75,
                        col='black', lwd=1, lty=1, main="Error / Dimension",
                        xlab="iteration", ylab="error / dimensionality", ...){
    bg.cols <- hsvScale(1:ncol(x$mapDims), alpha=bg.alpha, sat=bg.sat, val=bg.val)
    x.points <- 1:length(x$stress)
    max.x=length(x$stress)
    plot(x.points, x$stress, type='n', xlab=xlab,
         ylab=ylab, main=main, ...)
    y.range <- range(x$stress)
    y.span <- y.range[2] - y.range[1]
    max.dim <- sum(x$mapDims[1,])
    ## draw the background to indicate dimensionality
    for(i in ncol(x$mapDims):2){
        d <- apply(x$mapDims[, i:1], 1, sum)
        d <- y.range[1] + (d/max.dim)*y.span
        polygon( c(1, x.points, max.x), c(y.range[1], d, y.range[1]), col=bg.cols[i], border=NA )
        
    }
    d <- x$mapDims[,1]
    d <- y.range[1] + (d/max.dim) * y.span
    polygon( c(1, x.points, max.x), c(y.range[1], d, y.range[1]), col=bg.cols[1], border=NA )
    points(x.points, x$stress, type='l', lwd=lwd, col=col, lty=lty)
}

parallelDimFactors <- function(dim, iteration.no, red.end=iteration.no*0.75, target.dim=2){
    dimFactors <- matrix(nrow=iteration.no, ncol=dim, data=1.0)
    if(dim > target.dim){
        dimFactors[1:red.end, (target.dim+1):ncol(dimFactors) ] <- seq(from=1.0, to=0.0, length.out=red.end)
        dimFactors[(red.end+1):nrow(dimFactors), (target.dim+1):ncol(dimFactors) ] <- 0
    }
    dimFactors
}

parallelExpDimFactors <- function(dim, iteration.no, target.dim=2, red.end=iteration.no * 0.9){
    dimFactors <- matrix(nrow=iteration.no, ncol=dim, data=1.0)
    if(dim > target.dim){
        dimFactors[1:red.end, (target.dim+1):ncol(dimFactors) ] <- 2^( -seq(from=0, to=10, length.out=(red.end)) )
        dimFactors[(red.end+1):nrow(dimFactors), (target.dim+1):ncol(dimFactors) ] <- 0
    }
    dimFactors
}

serialDimFactors <- function(dim, iteration.no, red.end=iteration.no*0.75, target.dim=2){
    dimFactors <- matrix(nrow=iteration.no, ncol=dim, data=1.0)
    if(dim > target.dim){
        d.l <- as.integer(red.end / (dim - target.dim))
        red.i <- 1
        for(i in (dim):(target.dim+1)){
            dimFactors[red.i:(red.i+d.l-1), i] <- seq(from=1.0, to=0, length.out=d.l)
            dimFactors[(red.i+d.l):nrow(dimFactors), i] <- 0
            red.i <- red.i + d.l
        }
    }
    dimFactors
}

kStress <- function(sq){
    sqrt( sq$sq_stress/ sq$total_sq_dist )
}

nStress <- function(sq){
    sq$stress / sq$total_dist
}

