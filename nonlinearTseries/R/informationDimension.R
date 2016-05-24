################################################################################
#' Information dimension
#' @description
#' Functions for estimating the information dimension of a dynamical 
#' system from 1-dimensional time series using Takens' vectors
#' @details 
#' The information dimension is a particular case of the generalized correlation dimension
#' when setting the order q = 1. It is possible to demonstrate that the information dimension
#' \eqn{D_1}{D1} may be defined as:
#' \eqn{D_1=lim_{r \rightarrow 0} <\log p(r)>/\log(r)}{D1=lim{r->0} <ln p(r)>/ln(r)}.
#' Here, \eqn{p(r)} is the probability of finding a neighbour in a neighbourhood of size \eqn{r} and 
#' <> is the mean value. Thus, the information dimension specifies how the average
#' Shannon information scales with the radius \eqn{r}. The user should compute
#' the information dimension for different embedding dimensions for checking 
#' if \eqn{D_1}{D1} saturates.
#' 
#' In order to estimate \eqn{D_1}{D1}, the algorithm looks for the scaling behaviour of the the average
#' radius that contains a given portion (a "fixed-mass") of the total points in the phase space. By performing
#' a linear regression of \eqn{\log(p)\;Vs.\;\log(<r>)}{ln p Vs ln <r>} (being \eqn{p} the fixed-mass of the total points), an estimate
#' of \eqn{D_1}{D1} is obtained. 
#' 
#' The algorithm also introduces a variation of \eqn{p} for achieving a better performance: 
#' for small values of \eqn{p}, all the points in the time series (\eqn{N}) are considered for obtaining
#' \eqn{p=n/N}. Above a maximum number of neighbours \eqn{kMax}, the algorithm obtains \eqn{p} by decreasing the number
#' of points considerd  from the time series  \eqn{M<N}. Thus \eqn{p = kMax/M}.
#' 
#' Even with these improvements, the calculations for the information dimension are heavier than
#' those needed for the correlation dimension. 
#' @param time.series The original time series from which the information dimension will be estimated.
#' @param min.embedding.dim Integer denoting the minimum dimension in which we shall embed the time.series (see \link{buildTakens}). 
#' @param max.embedding.dim Integer denoting the maximum dimension in which we shall embed the time.series (see \link{buildTakens}).Thus,
#' we shall estimate the information dimension between \emph{min.embedding.dim} and \emph{max.embedding.dim}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors (see \code{\link{buildTakens}}).
#' @param min.fixed.mass Minimum percentage of the total points that the algorithm shall use for the estimation.
#' @param max.fixed.mass Maximum percentage of the total points that the algorithm shall use for the estimation.
#' @param number.fixed.mass.points The number of different \emph{fixed mass} fractions between \emph{min.fixed.mass}
#' and \emph{max.fixed.mass} that the algorithm will use for estimation.
#' @param radius Initial radius for searching neighbour points in the phase space. Ideally, it should be small
#' enough so that the fixed mass contained in this radius is slightly greater than the \emph{min.fixed.mass}. However,
#' whereas the radius is not too large (so that the performance decreases) the choice is not critical.
#' @param increasing.radius.factor Numeric value. If no enough neighbours are found within \emph{radius}, the radius
#' is increased by a factor \emph{increasing.radius.factor} until succesful. Default: sqrt(2) = 1.414214.
#' @param number.boxes Number of boxes that will be used in the box assisted algorithm (see \link{neighbourSearch}).
#' @param number.reference.vectors Number of reference points that the routine will try to use, saving computation time.
#' @param theiler.window Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  theiler.window time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations. 
#' @param kMax Maximum number of neighbours used for achieving p with all the points from the time series (see Details). 
#' @param do.plot Logical value. If TRUE (default value), a plot of the correlation sum is shown.
#' @param ... Additional graphical parameters.
#' @return A \emph{infDim} object that consist of a list with two components: \emph{log.radius} and \emph{fixed.mass}. \emph{log.radius} contains
#' the average log10(radius) in which the \emph{fixed.mass} can be found.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#'  @examples
#' \dontrun{
#' x=henon(n.sample=  3000,n.transient= 100, a = 1.4, b = 0.3, 
#'         start =  c(0.8253681, 0.6955566), do.plot = FALSE)$x
#' 
#' leps = infDim(x,min.embedding.dim=2,max.embedding.dim = 5,
#'               time.lag=1, min.fixed.mass=0.04, max.fixed.mass=0.2,
#'               number.fixed.mass.points=100, radius =0.001, 
#'               increasing.radius.factor = sqrt(2), number.boxes=100, 
#'               number.reference.vectors=100, theiler.window = 10, 
#'               kMax = 100,do.plot=FALSE)
#' 
#' plot(leps,type="l")
#' colors2=c("#999999", "#E69F00", "#56B4E9", "#009E73", 
#'           "#F0E442", "#0072B2", "#D55E00")
#' id.estimation = estimate(leps,do.plot=TRUE,use.embeddings = 3:5,
#'                          fit.lwd=2,fit.col=1,
#'                          col=colors2)
#' cat("Henon---> expected: 1.24    predicted: ", id.estimation ,"\n")
#' }
#' @rdname infDim
#' @export infDim
#' @exportClass infDim
#' @useDynLib nonlinearTseries
#' @seealso \code{\link{corrDim}}.
infDim <- 
  function(time.series, min.embedding.dim=2, 
           max.embedding.dim = min.embedding.dim,time.lag=1,
           min.fixed.mass, max.fixed.mass, number.fixed.mass.points = 10,
           radius, increasing.radius.factor = sqrt(2), number.boxes=NULL,
           number.reference.vectors=5000, theiler.window = 1,
           kMax = 1000,do.plot=TRUE,...){
   
    
    embeddings =  min.embedding.dim:max.embedding.dim
    n.embeddings = length(embeddings)
    fixed.mass.vector = 10^(seq(log10(min.fixed.mass),log10(max.fixed.mass),
                                length.out=number.fixed.mass.points))
    infDim.matrix = matrix(0, ncol = number.fixed.mass.points, nrow = n.embeddings)
    dimnames(infDim.matrix)  = list(embeddings,fixed.mass.vector)
    for (m in embeddings){
      infDim.result = 
        infDimSingleDimension(
          time.series, m, time.lag, min.fixed.mass,
          max.fixed.mass, number.fixed.mass.points, radius, 
          increasing.radius.factor, number.boxes,
          number.reference.vectors, theiler.window,
          kMax)
      infDim.matrix[as.character(m),] = infDim.result$lr
    }
    # Not using the correction of the log(p) suggested by Kantz,
    # the correction is stored in infDim.result$lfp.
    information.dimension.structure = list( fixed.mass = fixed.mass.vector, 
                                            log.radius =  infDim.matrix,
                                            embedding.dims = embeddings)
    class(information.dimension.structure) = "infDim"
    # add attributes
    id=deparse(substitute(time.series))
    attr(information.dimension.structure,"time.lag") = time.lag
    attr(information.dimension.structure,"id") = id
    attr(information.dimension.structure,"theiler.window") = theiler.window
    
    
    if (do.plot){
      plot(information.dimension.structure,...)
    }

    
    return(information.dimension.structure)
}


# Private function... computes the average log radius for a given fixed mass
# vector and a single embedding dimension.
infDimSingleDimension <- 
  function(time.series, embedding.dim, time.lag, min.fixed.mass,
                max.fixed.mass, number.fixed.mass.points, radius, 
                increasing.radius.factor, number.boxes,
                number.reference.vectors, theiler.window,
                kMax){
  #estimates number.boxes if it has not been specified
  takens = buildTakens(time.series,embedding.dim,time.lag)
  numberTakens = nrow(takens)
  embedding.dim = ncol(takens)
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  fixedMassVector = 10^(seq(log10(min.fixed.mass),log10(max.fixed.mass),length.out=number.fixed.mass.points))
  boxes=rep(0,number.boxes*number.boxes+1)
  averageLogRadiusVector = rep(0,number.fixed.mass.points)
  lnFixedMassVector = rep(0,number.fixed.mass.points)
  
  c.result=.C("informationDimension", takens = as.double(takens), 
              numberTakens = as.integer(numberTakens),
              embeddingD = as.integer(embedding.dim),
              fixedMassVector = as.double(fixedMassVector), 
              lnFixedMassVector = as.double(lnFixedMassVector), 
              fixedMassVectorLength = as.integer(number.fixed.mass.points),
              eps = as.double(radius),
              increasingEpsFactor = as.double(increasing.radius.factor),
              numberBoxes = as.integer(number.boxes), boxes = as.integer(boxes), 
              numberReferenceVectors = as.integer( number.reference.vectors), 
              theilerWindow = as.integer(theiler.window), kMax = as.integer(kMax),
              averageLogRadiusVector = as.double(averageLogRadiusVector),
              PACKAGE="nonlinearTseries")
 
  
  return( list(lr = c.result$averageLogRadiusVector,
              lfp = c.result$lnFixedMassVector))
}

#' Obtain the fixed mass vector used in the information dimension algorithm.
#' @param x A \emph{infDim} object.
#' @return A numeric vector representing the fixed mass vector used
#' in the information dimension algorithm represented by the \emph{infDim} object.
#' @seealso \code{\link{infDim}}
#' @export fixedMass
fixedMass = function(x){
  UseMethod("fixedMass")
}

#' @return The \emph{fixedMass} function returns the fixed mass vector used
#' in the information dimension algorithm.
#' @rdname infDim
#' @method fixedMass infDim
#' @export
fixedMass.infDim = function(x){
  return (x$fixed.mass)
}


#' Obtain the the average log(radius) computed
#' on the information dimension algorithm.
#' @param x A \emph{infDim} object.
#' @return A numeric vector representing the average log(radius) computed
#' on the information dimension algorithm represented by the \emph{infDim} object.
#' @seealso \code{\link{infDim}}
#' @export logRadius
logRadius = function(x){
  UseMethod("logRadius")
}


#' @return The \emph{logRadius} function returns the average log(radius) computed
#' on the information dimension algorithm.
#' @rdname infDim
#' @method fixedMass infDim
#' @export
logRadius.infDim = function(x){
  return (x$log.radius)
}


#' @return The \emph{embeddingDims} function returns the 
#' embeddings in which the information dimension was computed
#' @rdname infDim
#' @method embeddingDims infDim
#' @export
embeddingDims.infDim = function(x){
  return (embeddingDims.default(x))
}

#' @return The 'estimate' function estimates the information dimension of the 
#' 'infDim' object by by averaging the slopes of the
#' embedding dimensions specified in the \emph{use.embeddings} parameter. The
#' slopes are determined  by performing a linear regression
#' over the fixed mass' range specified in 'regression.range'. If do.plot is TRUE,
#' a graphic of the regression over the data is shown.
#' @param x A \emph{infDim} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @param use.embeddings A numeric vector specifying which embedding dimensions should 
#' the \emph{estimate} function use to compute the information dimension.
#' @param fit.col A vector of colors to plot the regression lines.
#' @param fit.lty The type of line to plot the regression lines.
#' @param fit.lwd The width of the line for the regression lines.
#' @param lty The line type of the <log10(radius)> functions.
#' @param lwd The line width of the <log10(radius)> functions.
#' @rdname infDim
#' @export
#' @method estimate infDim
#' @method estimate
estimate.infDim = function(x, regression.range=NULL, do.plot=TRUE,
                           use.embeddings = NULL,
                           col=NULL,pch=NULL,
                           fit.col=NULL, fit.lty=2,fit.lwd=2,
                           add.legend=T, lty=1,lwd=1,...){
  if (is.null(regression.range)){
    min.fixed.mass = min(x$fixed.mass) # the first position is always 0
    max.fixed.mass = max(x$fixed.mass)
  }else{
    min.fixed.mass = regression.range[[1]]
    max.fixed.mass = regression.range[[2]]
  }
  if (!is.null(use.embeddings)){
    log.radius = x$log.radius[as.character(use.embeddings),]
  }else{
    use.embeddings = x$embedding.dims
    log.radius=x$log.radius
  }
  
  n.embeddings = length(use.embeddings)
  indx = which(x$fixed.mass>=min.fixed.mass & x$fixed.mass<=max.fixed.mass)
  x.values = log10(x$fixed.mass[indx])
  information.dimension=c()
  if (do.plot){
    if (add.legend){
      current.par =  par(no.readonly = TRUE)
      on.exit(par(current.par))
      layout(rbind(1,2), heights=c(8,2))
    }
    # eliminate unnecessary embedding dimensions
    reduced.x = x
    reduced.x$log.radius = NULL
    reduced.x$log.radius = log.radius
    reduced.x$embedding.dims = NULL
    reduced.x$embedding.dims = use.embeddings
    # obtain vector of graphical parameters if not specified
    col = vectorizePar(col,n.embeddings)
    pch = vectorizePar(pch,n.embeddings)
    fit.col = vectorizePar(fit.col,n.embeddings,col)
    plot(reduced.x,col=col,pch=pch,lty=lty,lwd=lwd,add.legend=F,
         localScalingExp=F,...)
  }
  for(i in 1:n.embeddings){
    y.values = log.radius[as.character(use.embeddings[[i]]),indx]
    reg = lm(y.values~x.values)
    information.dimension = c(information.dimension,1/reg$coefficients[[2]])
    # plotting
    if (do.plot){
      lines(x$fixed.mass[indx],reg$fitted.values,
            col=fit.col[[i]],lwd=fit.lwd,lty=fit.lty,...)
    }  
  }
  if (add.legend && do.plot ){
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("center","groups",ncol=ceiling(n.embeddings/2), 
           bty="n", col=col,lty=rep(1,n.embeddings),pch=pch,
           lwd=rep(2.5,n.embeddings),
           legend=use.embeddings, title="Embedding dimension")
  }  
  
  information.dimension = mean(information.dimension)
  return(information.dimension)
  
}


#' @return The 'plot' function plots the computations performed for the 
#' information dimension estimate:  a graphic of <log10(radius)> Vs fixed mass.
#' Additionally, the inverse of the local scaling exponents can be plotted.
#' @param main A title for the plot.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param type Type of plot (see \code{\link[graphics]{plot}}).
#' @param log A character string which contains "x" if the x axis is to be 
#' logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx" if both
#'  axes are to be logarithmic.
#' @param ylim Numeric vector of length 2, giving the y coordinates range.
#' @param col Vector of colors for each of the dimensions of the plot.
#' @param pch Vector of symbols for each of the dimensions of the plot.
#' @param localScalingExp add a plot of the local information dimension 
#' scaling exponents?
#' @param add.legend add a legend to the plot?
#' @rdname infDim
#' @export
#' @method plot infDim
#' @method plot
plot.infDim = function(x, main="Information Dimension",
                       xlab="fixed mass (p)", ylab="<log10(radius)>",
                       type="b", log="x",ylim=NULL, col=NULL,pch=NULL,
                       localScalingExp=T,  add.legend=T,...){
  if ( add.legend || localScalingExp ){
    current.par =  par(no.readonly = TRUE)
    on.exit(par(current.par))
  }
  if (add.legend && localScalingExp){
    # 3 regions
    layout(rbind(1,2,3), heights=c(4,4,2))
  }else{
    if (add.legend){
      # add legend
      layout(rbind(1,2), heights=c(8,2))
    }
    if (localScalingExp){
      # add local slopes
      layout(rbind(1,2), heights=c(5,5))
    }
  }
  n.embeddings = length(x$embedding.dims)
  if (is.null(ylim)) ylim=range(x$log.radius)
  col = vectorizePar(col,n.embeddings)
  pch = vectorizePar(pch,n.embeddings)
  
  plot(x$fixed.mass,x$log.radius[as.character(x$embedding.dims[[1]]),],
       type=type,log=log, col=col[[1]], pch=pch[[1]],ylim=ylim, xlab=xlab,
       ylab=ylab,main=main,... )  
  if (n.embeddings >=2){
    for (i in 2:n.embeddings){
      lines(x$fixed.mass,x$log.radius[as.character(x$embedding.dims[[i]]),],
            type=type,col=col[[i]], pch=pch[[i]],...)
    }
  }
  #local slopes
  if (localScalingExp){
    plotLocalScalingExp(x,type=type,log=log,col=col,pch=pch,xlab=xlab,
                        add.legend=F,...)
  }
  ### add legend
  if (add.legend){
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("center","groups",ncol=ceiling(n.embeddings/2), 
           bty="n", col=col,lty=rep(1,n.embeddings),pch=pch,
           lwd=rep(2.5,n.embeddings),
           legend=x$embedding.dims, title="Embedding dimension")
  }
}



#' @return The \emph{plotLocalScalingExp} function plots the inverse of the 
#' local scaling exponentes of the information dimension 
#' (for reasons of numerical stability).
#' @rdname infDim
#' @method plotLocalScalingExp infDim
#' @export
#' @method plotLocalScalingExp
plotLocalScalingExp.infDim =  function(x,main="Local scaling exponents d1(p)",
                                       xlab="fixed mass p",ylab="1/d1(p)",
                                       type="b",log="x",ylim=NULL,col=NULL,pch=NULL,
                                       add.legend=T,...){
  
  n.embeddings = length(x$embedding.dims)
  if (length(x$fixed.mass) <= 1){
    stop("Cannot compute local scaling exponents (not enough points in the information dimension matrix)")
  }
  
  if (add.legend){
    current.par =  par(no.readonly = TRUE)
    on.exit(par(current.par))
    layout(rbind(1,2), heights=c(8,2))
  }  
  
  lfm = log10(x$fixed.mass)
  derivative = t(apply(x$log.radius,MARGIN=1,differentiate,h=lfm[[2]] - lfm[[1]]))
  fixed.mass.axis = differentiateAxis(x$fixed.mass)
  
  col = vectorizePar(col,n.embeddings)
  pch = vectorizePar(pch,n.embeddings)
  if (is.null(ylim)) ylim=range(derivative)
  plot(fixed.mass.axis,derivative[1,],type=type,log=log,col=col[[1]],
       pch=pch[[1]],ylim=ylim,xlab=xlab,ylab=ylab,main=main,...) 
  if (n.embeddings >=2){
    for (i in 2:n.embeddings){
      lines(fixed.mass.axis,derivative[i,],type=type,col=col[[i]],
            pch=pch[[i]],...)
    }
  }
  if (add.legend){
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("center","groups",ncol=ceiling(n.embeddings/2), 
           bty="n", col=col,lty=rep(1,n.embeddings),pch=pch,
           lwd=rep(2.5,n.embeddings),
           legend=x$embedding.dims, title="Embedding dimension")
  }
  
}
