################################################################################
#' Maximum lyapunov exponent
#' @description
#' Functions for estimating the maximal Lyapunov exponent of a dynamical system
#'  from 1-dimensional time series using Takens' vectors.
#' @details It is a well-known fact that close trajectories diverge exponentially fast in a chaotic system. The 
#' averaged exponent that determines the divergence rate is called the Lyapunov exponent (usually denoted with \eqn{\lambda}{lambda}). 
#' If \eqn{\delta(0)}{delta(0)} is the distance between two Takens' vectors in the embedding.dim-dimensional space, we expect that the distance
#' after a time \eqn{t} between the two trajectories arising from this two vectors fulfills:
#' \deqn{\delta (n) \sim \delta (0)\cdot exp(\lambda \cdot t)}{\delta (n) is.approximately \delta (0) exp(\lambda *t).}
#' The lyapunov exponent is estimated using the slope obtained by performing a linear regression of 
#' \eqn{S(t)=\lambda \cdot t \sim log(\delta (t)/\delta (0))}{S(t)=\lambda *t is.approximately log(\delta (t)/\delta (0))} 
#' on  \eqn{t}. \eqn{S(t)} will be estimated by averaging the divergence of several reference points.
#' 
#' The user should plot \eqn{S(t) Vs t} when looking for the maximal lyapunov exponent and, if for some temporal range
#' \eqn{S(t)} shows a linear behaviour, its slope is an estimate of the maximal Lyapunov exponent per unit of time. The estimate
#'  routine allows the user to get always an estimate of the maximal Lyapunov exponent, but the user must check that there is a linear region in the  
#' \eqn{S(t) Vs t}. If such a region does not exist, the estimation should be discarded. The computations should be performed
#' for several embedding dimensions in order to check that the Lyapunov exponent does not depend on the embedding dimension.
#' @param time.series The original time series from which the maximal Lyapunov exponent will be estimated
#' @param min.embedding.dim Integer denoting the minimum dimension in which we shall embed the time.series (see \link{buildTakens}). 
#' @param max.embedding.dim Integer denoting the maximum dimension in which we shall embed the time.series (see \link{buildTakens}).Thus,
#' we shall estimate the Lyapunov exponent between \emph{min.embedding.dim} and \emph{max.embedding.dim}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors (see \link{buildTakens}).
#' @param radius Maximum distance in which will look for nearby trajectories.
#' @param theiler.window Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  \emph{theiler.window} time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations. 
#' @param min.neighs Minimum number of neighbours that a Takens' vector must have to be considered
#' a reference point.
#' @param min.ref.points Number of reference points that the routine will try to use. The routine stops when it finds 
#' \emph{min.ref.points} reference points, saving computation time.
#' @param max.time.steps Integer denoting the number of time steps marking the end of the linear region.
#' @param number.boxes Number of boxes that will be used in the box assisted algorithm (see \link{neighbourSearch}).
#' @param sampling.period Sampling period of the time series. When dealing with a discrete
#' system, the \emph{sampling.period} should be set to 1.
#' @param do.plot Logical value. If TRUE (default value), a plot of \eqn{S(t)} Vs  \eqn{t} is shown.
#' @param ... Additional plotting parameters. 
#'  @return A list with three components named  \eqn{time} and \eqn{s.function}.
#' \eqn{time} is a vector containing the temporal interval where the system 
#' evolves. It ranges from 0 to 
#' \emph{\eqn{max.time.steps \cdot sampling.period}{max.time.steps * sampling.period}}.
#' \eqn{s.function} is a matrix containing the 
#' values of the \eqn{S(t)} for each t in the time vector(the columns) and each 
#' embedding dimension (the rows).
#' @references  
#' Eckmann, Jean-Pierre and Kamphorst, S Oliffson and Ruelle, David and Ciliberto, S and others. Liapunov exponents from time series.
#' Physical Review A, 34-6, 4971--4979, (1986).
#' 
#' Rosenstein, Michael T and Collins, James J and De Luca, Carlo J.A practical method for calculating largest Lyapunov exponents from small data sets.
#' Physica D: Nonlinear Phenomena, 65-1, 117--134, (1993).
#' @examples \dontrun{
#' ## Henon System
#' h=henon(n.sample=  5000,n.transient= 100, a = 1.4, b = 0.3, 
#'         start = c(0.63954883, 0.04772637), do.plot = FALSE) 
#' my.ts=h$x 
#' ml=maxLyapunov(time.series=my.ts,
#'                min.embedding.dim=2,
#'                max.embedding.dim=5,
#'                time.lag=1,
#'                radius=0.001,theiler.window=4,
#'                min.neighs=2,min.ref.points=500,
#'                max.time.steps=40,do.plot=FALSE)
#' plot(ml)
#' ml.estimation = estimate(ml,regression.range = c(0,15),
#'                          use.embeddings=4:5,
#'                          do.plot = TRUE)
#' # The max Lyapunov exponent of the Henon system is 0.41
#' cat("expected: ",0.41," calculated: ",ml.estimation,"\n")
#' 
#' ## Rossler system
#' r=rossler(a=0.15,b=0.2,w=10,start=c(0,0,0), time=seq(0,1000,0.1),do.plot=FALSE)
#' my.ts=r$x
#' use.cols = c("#999999","#E69F00","#56B4E9")
#' ml=maxLyapunov(time.series=my.ts,min.embedding.dim=5,max.embedding.dim = 7,
#'                time.lag=12,radius=0.1,theiler.window=50,
#'                min.neighs=5,min.ref.points=length(r),
#'                max.time.steps=300,number.boxes=NULL,
#'                sampling.period=0.1,do.plot=TRUE,
#'                col=use.cols)
#' #  The max Lyapunov exponent of the Rossler system is 0.09
#' ml.est=estimate(ml,col=use.cols,do.plot=T,
#'                 fit.lty=1,
#'                 fit.lwd=5)
#' cat("expected: ",0.090," calculated: ",ml.est,"\n")
#' }
#' @author Constantino A. Garcia
#' @rdname maxLyapunov
#' @export maxLyapunov
#' @exportClass maxLyapunov
#' @useDynLib nonlinearTseries
maxLyapunov=function(time.series,
                     min.embedding.dim=2,
                     max.embedding.dim=min.embedding.dim,
                     time.lag=1,radius,theiler.window=1,min.neighs=5,
                     min.ref.points=500,max.time.steps=10,number.boxes=NULL,
                     sampling.period=1,do.plot=TRUE,...){
  
  embeddings = min.embedding.dim:max.embedding.dim
  n.embeddings = length(embeddings)
  s.function = matrix(0, ncol = (max.time.steps+1), nrow = n.embeddings)
  dimnames(s.function) = list(embeddings,0:max.time.steps)
  for (m in embeddings){
    s.function[as.character(m),] =
      maxLyapunovSingleDimension(time.series,takens = NULL,
      m, time.lag, radius,theiler.window,
      min.neighs,min.ref.points,
      max.time.steps,number.boxes,
      sampling.period)
  }
  
  time=(0:max.time.steps)*sampling.period
  max.lyapunov.structure = list(time=time,s.function=s.function, 
                                embedding.dims = min.embedding.dim:max.embedding.dim)
  
  class(max.lyapunov.structure) = "maxLyapunov"
  id=deparse(substitute(time.series))
  attr(max.lyapunov.structure,"id") = id
  attr(max.lyapunov.structure,"time.lag") = time.lag
  attr(max.lyapunov.structure,"theiler.window") = theiler.window  
  attr(max.lyapunov.structure,"min.neighs") = min.neighs
  attr(max.lyapunov.structure,"min.ref.points") = min.ref.points
  
  #plot
  if (do.plot){
    plot(max.lyapunov.structure,...)
  }
  #return results
  return (max.lyapunov.structure)
  
}

# Private function for the computation of the s.function of a given 
# embedding dimension
maxLyapunovSingleDimension=function(time.series,takens,embedding.dim,time.lag,
                                    radius,theiler.window,min.neighs,
                                    min.ref.points,max.time.steps,number.boxes,
                                    sampling.period,do.plot){
  #C parameters
  if (is.null(takens)) takens=buildTakens(time.series,embedding.dim=embedding.dim,time.lag=time.lag)
  if (is.null(number.boxes)) {
      number.boxes = estimateNumberBoxes(time.series, radius)
  }
  numberTakens=nrow(takens)
  Sdn=rep(0,max.time.steps+1)
  
  sol=.C("maxLyapunov",timeSeries=as.double(time.series),takens=as.double(takens),
         tau=as.integer(time.lag),numberTakens=as.integer(numberTakens),
         embeddingD=as.integer(embedding.dim),eps=as.double(radius),
         Sdn=as.double(Sdn),nmax=as.integer(max.time.steps),
         nminRP=as.integer(min.ref.points),neighMin=as.integer(min.neighs),
         numberBoxes=as.integer(number.boxes),tdist=as.integer(theiler.window),
         PACKAGE="nonlinearTseries")
  
  # create lyapunov structure
  return(sol$Sdn)
}

#' Returns the time in which the divergence of close trajectories was computed 
#' in order to estimate the maximum Lyapunov exponent.
#' @param x A \emph{maxLyapunov} object.
#' @return A numeric vector representing the time in which the divergence of
#' close trajectories was computed.
#' @seealso \code{\link{maxLyapunov}}
#' @export divTime
divTime = function(x){
  UseMethod("divTime")
}


#' @return The \emph{divTime} function returns the 
#' time in which the divergence of close trajectories was computed.
#' @rdname maxLyapunov
#' @method divTime maxLyapunov
#' @export
divTime.maxLyapunov = function(x){
  return (x$time)
}

#' @return The \emph{embeddingDims} function returns the 
#' embeddings in which the divergence of close trajectories was computed
#' @rdname maxLyapunov
#' @method embeddingDims maxLyapunov
#' @export
embeddingDims.maxLyapunov = function(x){
  return (embeddingDims.default(x))
}

#' Returns the rate of divergence of close trajectories needed for the maximum Lyapunov
#' exponent estimation.
#' @param x A \emph{maxLyapunov} object.
#' @return A numeric matrix representing the time in which the divergence of
#' close trajectories was computed. Each row represents an embedding dimension
#' whereas that each column represents an specific moment of time.
#' @seealso \code{\link{maxLyapunov}}
#' @export divergence
divergence = function(x){
  UseMethod("divergence")
}


#' @return The \emph{divergence} function returns the 
#' rate of divergence of close trajectories needed for the maximum Lyapunov
#' exponent estimation.
#' @rdname maxLyapunov
#' @method divergence maxLyapunov
#' @export
divergence.maxLyapunov = function(x){
  return (x$s.function)
}

#' @rdname maxLyapunov
#' @method plot maxLyapunov  
#' @param main A title for the plot.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param type Type of plot (see \code{\link[graphics]{plot}}).
#' @param col Vector of colors for each of the dimensions of the plot.
#' @param pch Vector of symbols for each of the dimensions of the plot.
#' @param add.legend add a legend to the plot?
#' @method plot maxLyapunov
#' @export
#' @method plot
plot.maxLyapunov= function(x, main="Estimating maximal Lyapunov exponent",
                           xlab="time t",ylab="S(t)",type="p",col=NULL,
                           pch=NULL, add.legend=T,...){
  embeddings = x$embedding.dims
  n.embeddings = length(embeddings)
  # obtain default parameter vectors if not specified
  col = vectorizePar(col,n.embeddings)
  pch = vectorizePar(pch,n.embeddings)
  for (i in 1:n.embeddings){
    if (i!=1){
      lines(x$time[-1],x$s.function[as.character(embeddings[[i]]),-1],
            type=type, col = col[[i]],pch=pch[[i]],...)
    }else{# first time: use plot
      plot(x$time[-1],x$s.function[as.character(embeddings[[i]]),-1],
           xlab=xlab,ylab=ylab, main=main,type=type, col=col[[1]],
           pch=pch[[1]],...)    
    }
  }
  if (add.legend){
    legend("bottomright",
           col=col[1:n.embeddings],
           pch=pch[1:n.embeddings],lty=rep(1,n.embeddings), 
           lwd=rep(2.5,n.embeddings),bty="n",
           legend=embeddings,title="Embedding dimension")
  }
  
}

#' @return In order to obtain an estimation of the Lyapunov exponent the user can use the
#' \emph{estimate} function. The  \eqn{estimate} function allows the user to obtain
#' an estimation of the maximal Lyapunov exponent by averaging the slopes of the
#' embedding dimensions specified in the \emph{use.embeddings} parameter. The
#' slopes are determined by performing a linear regression
#' over the radius' range specified in \emph{regression.range}
#' @param ylim Numeric vector of length 2, giving the y coordinates range.
#' @param x A \emph{maxLyapunov} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @param use.embeddings A numeric vector specifying which embedding dimensions should 
#' the \emph{estimate} function use to compute the Lyapunov exponent.
#' @param fit.col A vector of colors to plot the regression lines.
#' @param fit.lty The type of line to plot the regression lines.
#' @param fit.lwd The width of the line for the regression lines.
#' @rdname maxLyapunov
#' @export
#' @method estimate maxLyapunov
#' @method estimate
estimate.maxLyapunov= function(x,regression.range = NULL,
                               do.plot=FALSE,use.embeddings = NULL,
                               main="Estimating maximal Lyapunov exponent",
                               xlab="time t",ylab="S(t)",type="p",col=NULL,
                               pch=NULL, ylim=NULL, fit.col=NULL,fit.lty=2,
                               fit.lwd=2, add.legend=T,...){
  if (is.null(regression.range)){
    min.time = x$time[[2]] # the first position is always 0
    max.time = tail(x$time,1)
  }else{
    min.time = regression.range[[1]]
    max.time = regression.range[[2]]
  }
  if (!is.null(use.embeddings)){
    s.function = x$s.function[as.character(use.embeddings),]
  }else{
    use.embeddings = x$embedding.dims
    s.function=x$s.function
  }
  
  indx = which(x$time >= min.time & x$time <= max.time )
  x.values = x$time[indx]
  lyapunov.estimate = c()
  n.embeddings = length(use.embeddings)
  # color and symbols for the plots if needed
  if (do.plot){
    # obtain vector of graphical parameters if not specified
    col = vectorizePar(col,n.embeddings)
    pch = vectorizePar(pch,n.embeddings)
    fit.col = vectorizePar(fit.col,n.embeddings,col)  
  }
  for (i in 1:n.embeddings){
    current.embedding = use.embeddings[[i]]
    y.values = s.function[as.character(current.embedding),indx]
    fit=lm(y.values~x.values)
    lyapunov.estimate=c(lyapunov.estimate,fit$coefficients[[2]])
    if (do.plot){
      if (i!=1){
        lines(x$time,s.function[as.character(current.embedding),], 
              type=type, col = col[[i]],pch=pch[[i]],...)
      }else{
        if (is.null(ylim)) ylim=range(s.function)
        plot(x$time,s.function[as.character(current.embedding),],
             xlab=xlab,ylab=ylab, main=main,type=type, col=col[[1]],
             pch=pch[[1]], ylim = ylim,...)  
      }
      lines(x.values,fit$fitted.values,col=fit.col[[i]],
            lty=fit.lty,lwd=fit.lwd)
    }
  }
  lyapunov.estimate = mean(lyapunov.estimate)
  if (do.plot && add.legend){
    legend("bottomright",col=col[1:n.embeddings],
           pch=pch[1:n.embeddings],
           lty=rep(1,n.embeddings),bty="n", 
           lwd=rep(2.5,n.embeddings), legend=use.embeddings,
           title="Embedding dimension")
  }
  
  
  return (lyapunov.estimate)
}

