################################################################################
#' Sample Entropy (also known as Kolgomorov-Sinai Entropy)
#' @description
#' The Sample Entropy measures the complexity of a time series. Large values of 
#' the Sample Entropy indicate high complexity whereas that smaller values characterize
#' more regular signals.
#' @details  The sample entropy is computed using:
#' \deqn{h_q(m,r) = log(C_q(m,r)/C_{q}(m+1,r))}{hq(m,r) = log(Cq(m,r)/Cq(m+1,r)),}
#' where \emph{m} is the embedding dimension and \emph{r} is the radius of the neighbourhood. When 
#' computing the correlation dimensions we use the linear regions from the correlation
#' sums in order to do the estimates. Similarly, the sample entropy \eqn{h_q(m,r)}{hq(m,r)} 
#' should not change for both various \emph{m} and \emph{r}.
#' @param corrDim.object A \emph{corrDim} object from which the Sample Entropy
#' of the time series characterized by \emph{corrDim} shall be estimated.
#' @param do.plot do.plot Logical value. If TRUE (default value), a plot of the sample entropy is shown.
#' @param ... Additional plotting arguments.
#' @return A \emph{sampleEntropy} object that contains a list storing the sample entropy (\emph{sample.entropy}),
#' the embedding dimensions ( \emph{embedding.dims}) and radius (\emph{radius}) for which the sample entropy has 
#' been computed, and the order of the sample entropy (\emph{entr.order}). The sample entropy
#' is stored as a matrix in which each row contains the computations for a given embedding dimension and 
#' each column stores the computations for a given radius.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @examples
#' \dontrun{
#' x=henon(n.sample = 15000, n.transient = 100, a = 1.4, b = 0.3, 
#'         start = c(0.78,0.8165), do.plot = FALSE)$x
#' 
#' cd=corrDim(time.series=x,
#'            min.embedding.dim=2,max.embedding.dim=9,
#'            corr.order=2,time.lag=1,
#'            min.radius=0.05,max.radius=1,
#'            n.points.radius=100,
#'            theiler.window=20,
#'            do.plot=TRUE)
#' 
#' use.col = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
#'             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#' se=sampleEntropy(cd,do.plot=TRUE,col=use.col,
#'                  type="l",xlim=c(0.1,1),
#'                  add.legend=T)
#' se.est = estimate(se,
#'                   regression.range = c(0.4,0.6),
#'                   use.embeddings = 6:9,col=use.col,type="b")
#' print(se.est)
#' cat("Expected K2 = ",0.325," Estimated = ",mean(se.est),"\n")
#' }
#' @author Constantino A. Garcia
#' @rdname sampleEntropy
#' @export sampleEntropy
#' @exportClass sampleEntropy
sampleEntropy = function (corrDim.object, do.plot=TRUE,...){ 
  if (!inherits(corrDim.object, "corrDim")){
    stop("corrDim.object should be of class corrDim")
  }
  radius = radius(corrDim.object)
  corr.matrix = corrMatrix(corrDim.object)
  embeddings = embeddingDims(corrDim.object)
  number.embeddings = length(embeddings) - 1
  entropy = matrix(0,nrow= number.embeddings,ncol=length(radius))
  for (i in 1:number.embeddings){
    entropy[i,] = log(corr.matrix[i,]/corr.matrix[i+1,])
  }
  dimnames(entropy)=list(head(embeddings,-1),radius)
  sample.entropy = list(sample.entropy = entropy,embedding.dims = head(embeddings,-1),
                        entr.order=nlOrder(corrDim.object), radius=radius)
  class(sample.entropy)="sampleEntropy"
  
  attr(sample.entropy,"time.lag") = attr(corrDim.object,"time.lag") 
  attr(sample.entropy,"id") = attr(corrDim.object,"id") 
  attr(sample.entropy,"theiler.window") = attr(corrDim.object,"theiler.window") 
  
  if (do.plot){
   plot(sample.entropy,...) 
  }
  
  return (sample.entropy)

}

#' Returns the sample entropy function \eqn{h_q(m,r)} used for the computations
#' of the sample entropy.
#' @param x A \emph{sampleEntropy} object.
#' @return A numeric matrix representing the sample entropy function
#' \eqn{h_q(m,r)} obtained in #' the sample entropy computations represented
#' by the \emph{sampleEntropy} object.
#' @seealso \code{\link{sampleEntropy}}
#' @export sampleEntropyFunction
sampleEntropyFunction = function(x){
  UseMethod("sampleEntropyFunction")
}

#' @return The \emph{sampleEntropyFunction} returns the sample entropy function
#' \eqn{h_q(m,r)} used for the computations. The sample
#' entropy function is represented by a matrix. Each row represents a given
#' embedding dimension whereas that each column representes a different radius.
#' @rdname sampleEntropy
#' @method sampleEntropyFunction sampleEntropy
#' @export
sampleEntropyFunction.sampleEntropy = function(x){
  return (x$sample.entropy)
}

#' @return The \emph{nlOrder} function returns the order of the sample entropy.
#' @rdname sampleEntropy
#' @method nlOrder sampleEntropy
#' @export
nlOrder.sampleEntropy = function(x){
  return(x$entr.order)
}

#' @return The \emph{radius} function returns the radius on which the sample entropy
#'  function has been evaluated.
#' @rdname sampleEntropy
#' @method radius sampleEntropy
#' @export
radius.sampleEntropy = function(x){
  return (radius.default(x))
}

#' @return The \emph{embeddingDims} function returns the embedding dimensions 
#' on which the sample entropy function has been evaluated.
#' @rdname sampleEntropy
#' @method embeddingDims sampleEntropy
#' @export
embeddingDims.sampleEntropy = function(x){
  return (embeddingDims.default(x))
}


#' @return The \emph{plot} function shows the graphics for the sample entropy.
#' @rdname sampleEntropy
#' @method plot sampleEntropy
#' @param main A title for the plot.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param type Type of plot (see \code{\link[graphics]{plot}}).
#' @param col Vector of colors for each of the dimensions of the plot.
#' @param pch Vector of symbols for each of the dimensions of the plot
#' @param ylim Numeric vector of length 2, giving the y coordinates range..
#' @param add.legend add a legend to the plot?
#' @export
plot.sampleEntropy = function(x,main=NULL,xlab=NULL,ylab=NULL,type="l",
                              col=NULL, pch=NULL, ylim=NULL,add.legend=T,...){
  if (is.null(xlab)) xlab = expression("ln("*epsilon*")")
  if (is.null(ylab)) ylab = expression(h[q]*"("*epsilon*")")
  if (is.null(main)) main = 
    bquote("Sample entropy (q = "*.(x$entr.order)*")    "*h[.(x$entr.order)]*"("*epsilon*")")
  number.embeddings = length(x$embedding.dims)
  
  # obtain vector of graphical parameters if not specified
  col = vectorizePar(col,number.embeddings)
  pch = vectorizePar(pch,number.embeddings)
 
  if (add.legend){
    current.par =  par(no.readonly = TRUE)
    on.exit(par(current.par))
    layout(rbind(1,2), heights=c(8,2))
  }
  for (i in 1:number.embeddings){
     if (i == 1) {
       if (is.null(ylim)) ylim=range(x$sample.entropy)
       plot(x$radius,x$sample.entropy[1,],
            xlab = xlab,ylab = ylab,main=main,type=type,pch=pch[[i]],
            col=col[[i]],ylim=ylim,...)
     }else{
       lines(x$radius,x$sample.entropy[i,],type=type,pch=pch[[i]],
             col=col[[i]],...)
     }
  }
  if(add.legend){
    par(mar=c(0, 0, 0, 0))
    # c(bottom, left, top, right)
    plot.new()
    legend("center","groups",ncol=ceiling(number.embeddings/2),
           col=col,pch=pch,lty=rep(1,number.embeddings),
           lwd=rep(2.5,number.embeddings),bty="n",
           legend=x$embedding.dims,title="Embedding dimension")
  }
  
}


#' @details For each embedding dimension the sample
#' entropy is estimated by averaging  \deqn{h_q(m,r) = log(C_q(m,r)/C_{q}(m+1,r))}{hq(m,r) = log(Cq(m,r)/Cq(m+1,r))}
#' over the range specified by \emph{regression range} in the \emph{estimate} function.
#' @param x A \emph{sampleEntropy} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @param use.embeddings A numeric vector specifying which embedding dimensions should the \emph{estimate} function use to compute
#' the sample entropy.
#' @param fit.lty The type of line to plot the regression lines. 
#' @param fit.lwd	The width of the line for the regression lines.
#' @param fit.col A vector of colors to plot the regression lines.
#' @return The  \emph{estimate} function returns a vector storing the sample entropy estimate for each embedding dimension.
#' @rdname sampleEntropy
#' @method estimate sampleEntropy
#' @export
estimate.sampleEntropy = function(x,regression.range=NULL,do.plot=TRUE,
                                  use.embeddings=NULL,fit.col=NULL,
                                  fit.lty=2, fit.lwd=2,add.legend=T,...){
  if (is.null(regression.range)){
    regression.range = range(x$radius)
  }  
  if (is.null(use.embeddings)){
    use.embeddings = x$embedding.dims
  }
  #select only embeddings used in the object
  use.embeddings = intersect(as.numeric(use.embeddings),x$embedding.dims)
  len.use.embeddings = length(use.embeddings)
  if(len.use.embeddings<1) {stop("The embeddings specified are not present!\n")}
  #range
  indx = which(x$radius >= regression.range[[1]] & x$radius <=regression.range[[2]])
    
  if (len.use.embeddings == 1){
    sample.entropy.estimate = mean(x$sample.entropy[as.character(use.embeddings),indx])
  }else{
    sample.entropy.estimate = apply(x$sample.entropy[as.character(use.embeddings),indx],MARGIN=1,FUN=mean)  
  }  
  if (do.plot){
    #create new sample entropy object with the embedding.dims and radius range used
    x = list(sample.entropy = x$sample.entropy[as.character(use.embeddings),],
             embedding.dims = use.embeddings,entr.order=x$entr.order, radius=x$radius)
    class(x)="sampleEntropy"
    plotSampleEntropyEstimate(x,sample.entropy.estimate,
                              fit.col=fit.col,fit.lty=fit.lty,
                              fit.lwd=fit.lwd,add.legend=add.legend,...)
  }
  return(sample.entropy.estimate)
}

plotSampleEntropyEstimate = function(sampleEntropy.object,
                                     sample.entropy.estimate,main=NULL,
                                     xlab=NULL,ylab=NULL,type="l",pch=NULL,
                                     ylim=NULL,col=NULL,
                                     fit.col,fit.lty,fit.lwd,
                                     add.legend=T,...){
  
  number.embeddings = length(sampleEntropy.object$embedding.dims)
  
  # obtain vector of graphical parameters if not specified
  col = vectorizePar(col,number.embeddings)
  pch = vectorizePar(pch,number.embeddings)
  fit.col = vectorizePar(fit.col,number.embeddings,col)
  
  if (add.legend){
    current.par =  par(no.readonly = TRUE)
    on.exit(par(current.par))
    layout(rbind(1,2), heights=c(8,2))
  }
  plot(sampleEntropy.object,main=main,
       xlab=xlab,ylab=ylab,type=type,pch=pch,
       ylim=ylim,col=col,add.legend=F,...)
  for (i in 1:number.embeddings){
    if (i == 1) {
      abline(h=sample.entropy.estimate[[i]],
             col=fit.col[[i]],lty=fit.lty,lwd=fit.lwd)
    }else{
      abline(h=sample.entropy.estimate[[i]],
             col=fit.col[[i]],lty=fit.lty,lwd=fit.lwd)
    }
  }
  if(add.legend){
    par(mar=c(0, 0, 0, 0))
    # c(bottom, left, top, right)
    plot.new()
    legend("center","groups",ncol=ceiling(number.embeddings/2),
           col=col,pch=pch,lty=rep(1,number.embeddings),
           lwd=rep(2.5,number.embeddings),bty="n",
           legend=sampleEntropy.object$embedding.dims,
           title="Embedding dimension")
  }
  
}
