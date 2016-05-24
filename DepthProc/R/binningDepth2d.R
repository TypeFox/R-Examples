#'@title 2d Binning
#'@importFrom sm binning
#'@description A robust method of decreasing a sample size and therefore a complexity of a statistical procedure. The method may be used within a kernel density or a predictive distribution estimation.
#'
#'  @param x bivariate matrix containing data. Each row is viewed as one two-dimensional observation.
#'  @param binmethod A method for calculation center and dispersion measures. "LocDepth" uses location-scale depth, ?MAD? uses median and MAD in each dimension.
#'  @param nbins number of bins in each dimension
#'  @param k responsible for tightness of bins.
#'  @param remove_borders Logical, include or not marginal bins
#'  @param ... other arguments passed to depthMedian
#'  
#'  @return freq: a matrix containing the binned frequencies
#'  @return mid_x: mid points for x
#'  @return mid_y: mid points for y
#'  @return breaks_x: breaks for x
#'  @return breaks_y: breaks for y
#'  @return input_data: max_depth_x and  max_depth_y:
#'    
#'  
#'  
#' 
#'  @seealso \code{\link{depth}}
#'  @export
#'  
#'  @details
#'  
#'  Let us recall, that binning is a popular method of decreasing a sample size. To bin a window of  \eqn{ n }  points  \eqn{ {W}_{i,n}=\left\{{X}_{i-n+1},...,{X}_{i} \right\} }  to a grid  \eqn{ {{{X}'}_{1}},...,{{{X}'}_{m}} }  we simply assign each sample point  \eqn{ {{X}_{i}} }  to the nearest grid point  \eqn{ {{{X}'}_{j}} } . When binning is completed, each grid point  \eqn{ {{X}'}_{j} }  has an associated number  \eqn{ {c}_{i} } , which is the sum of all the points that have been assigned to  \eqn{ {{X}'}_{j} } . This procedure replaces the data  \eqn{ {W}_{i,n}=\left\{ {X}_{i-n+1},...,{X}_{i} \right\} }  with the smaller set  \eqn{ {{W}'}_{j,m}=\left\{ {{X}'}_{j-m+1},...,{{X}'}_{j} \right\} } . Although simple binning can speed up the computation, it is criticized for a lack of a precise approximate control over the accuracy of the approximation. Robust binning however stresses properties of the majority of the data and decreases the computational complexity of the DSA at the same time.                                                                                                                                                                                                                                                                                     
#'  
#'  For a 1D window  \eqn{ {W}_{i,n} } , let  \eqn{ {Z}_{i,n-k} }  denote a 2D window created basing on  \eqn{ {W}_{i,n} }  and consisted of  \eqn{ n-k }  pairs of observations and the  \eqn{ k }  lagged observations  \eqn{ {Z}_{i,n-k} } = \eqn{ \left\{ ({X}_{i-n-k},{X}_{i-n+1})\right\} } ,  \eqn{ 1\le i\le n-k }  . Robust 2D binning of the  \eqn{ {Z}_{i,n-p} }  is a very useful technique in a context of robust estimation of the predictive distribution of a time series (see \cite{Kosiorowski:2013b}).                                                                                                                                                                                                                                                                                     
#'  
#'  Assume we analyze a data stream  \eqn{ \{{X}_{t}\} }  using a moving window of a fixed length  \eqn{ n } , i.e.,  \eqn{ {W}_{i,n} }  and the derivative window  \eqn{ {Z}_{i,n-1} } . In a first step we calculate the weighted sample  \eqn{ L^p }  depth for  \eqn{ {W}_{i,n} } . Next we choose equally spaced grid of points  \eqn{ {l}_{1},...,{l}_{m} }  in this way that  \eqn{ [{{l}_{1}},{{l}_{m}}]\times [{{l}_{1}},{{l}_{m}}] }  covers fraction of the   \eqn{ \beta }  central points of  \eqn{ {Z}_{i,n-1} }  w.r.t. the calculated  \eqn{ L^p }  depth, i.e., it covers  \eqn{ {R}^{\beta }({Z}_{i,n-1}) }  for certain prefixed threshold  \eqn{ \beta \in (0,1) } . For both  \eqn{ {X}_{t} }  and  \eqn{ {X}_{t-1} }  we perform a simple binning using following bins:  \eqn{ (-\infty ,{l}_{1}) } ,  \eqn{ ({l}_{1},{l}_{2}) } ,...,  \eqn{ ({l}_{m},\infty ) } .                                                                                                                                                                                                                                                                                     
#'  
#'  For robust binning we reject "border" classes and further use only midpoints and binned frequencies for classes  \eqn{ ({l}_{1},{l}_{2}) } ,  \eqn{ ({l}_{2},{l}_{3}) } ,...,   \eqn{ ({l}_{m-1},{l}_{m}) } .
#'  
#'  @author Daniel Kosiorowski and Zygmunt Zawadzki from Cracow University of Economics.
#'  
#'  @references
#'  
#'  Hall, P., Wand, M. P. (1996) On the Accuracy of Binned Kernel Density Estimators, Journal of Multivariate Analysis archive, Volume 56 Issue 2, 165 - 184
#'  
#' Holmstrom, L. (2000) The Accuracy and the Computational Complexity of a Multivariate Binned Kernel Density Estimator, Journal of Multivariate Analysis, Volume 72, Issue 2, 264-309, http://dx.doi.org/10.1006/jmva.1999.1863. (http://www.sciencedirect.com/science/article/pii/S0047259X99918638)
#'  
#'  @examples
#'  
#' #EXAMPLE 1
#' Sigma1 = matrix(c(10,3,3,2),2,2)
#' X1 = mvrnorm(n= 8500, mu= c(0,0),Sigma1)
#' Sigma2 = matrix(c(10,0,0,2),2,2)
#' X2 = mvrnorm(n= 1500, mu= c(-10,6),Sigma2)
#' BALLOT = rbind(X1,X2)
#' train = sample(1:10000, 500)
#' data =BALLOT[train,]
#' plot(data)
#' 
#' b1=binningDepth2D(data, remove_borders = FALSE, nbins = 12, k = 1 )
#' b2=binningDepth2D(data, nbins = 12, k = 1,remove_borders = TRUE )
#' plot(b1)
#' plot(b2)
#' 
#' #EXAMPLE 2
#' data(under5.mort)
#' data(maesles.imm)
#' data2011=cbind(under5.mort[,22],maesles.imm[,22])
#' plot(binningDepth2D(data2011, nbins = 8, k = 0.5,remove_borders = TRUE ))
#'  
#'  @keywords
#'  multivariate
#'  nonparametric
#'  robust
#'  depth function
binningDepth2D = function(x, binmethod = "LocDepth", nbins = 8, k = 1, remove_borders = FALSE,  ...)
{
  createBin = function(x, nbins, mean = NULL)
  {
    if(binmethod == "LocDepth")
    {
      #if(!devel) dep_stat = sample.max.depth(as.numeric(x)) 
      #else 
      #{
        dep_stat = lsdSampleMaxDepth(x)
      #}
      mean = dep_stat@mu#["mu"]
      sigma = k*dep_stat@sigma#["sigma"]
      dep_stat = c(dep_stat@max_depth,dep_stat@mu, dep_stat@sigma)    
    }
    if(binmethod == "LP")
    {
      sigma  = mad(x)
      dep_stat = c(mu = mean, sigma = sigma)
    }
    
    range = range(x)
    
    if(nbins == "auto")
    {
      n_upper = 1:ceiling(abs(range[2]-mean)/sigma)*sigma+mean
      n_lower = -(ceiling(abs(range[1]-mean)/sigma):1)*sigma+mean
      breaks = c(n_lower,mean, n_upper)
    } else
    {
      #if(!remove_borders) 
      nbins = nbins-2
      
      if(nbins>0)
      {
        s_bound = 1:(nbins/2)*sigma
      }
      0:nbins*sigma
      breaks = sort(c(-Inf,-s_bound+mean,mean,s_bound+mean,Inf))
    }
    
    cut = cut(x, breaks=breaks)
    midpoints = sapply(2:length(breaks), function(x) mean(breaks[(x-1):x]))
    
    if(nbins != "auto")
    {
      midpoints[1] =  midpoints[2]-2*sigma #mean(x[x<breaks[2]])
      midpoints[length(midpoints)] =  midpoints[(length(midpoints)-1)]+2*sigma #mean(x[x>breaks[length(breaks)-1]])
      if(is.na(midpoints[1])) midpoints[1] = midpoints[2]-2*sigma
      if(is.na(midpoints[length(midpoints)])) midpoints[(length(midpoints))] = midpoints[(length(midpoints)-1)]+2*sigma
    }
    
    res = list(breaks, midpoints, dep_stat)
    return(res)
  }
  
  means = c(0,0)
  if(binmethod == "LP") means = depthMedian(x,method = "LP",...)
  
  tmp1 = createBin(x[,1], nbins, means[1])
  tmp2 = createBin(x[,2], nbins, means[2])
  
  
  b = cbind(tmp1[[1]],tmp2[[1]])
  b[b == Inf] = 1e6
  b[b == -Inf] = -1e6
  tmp = binning(x=x,breaks=b)$table.freq
  
  if(remove_borders == TRUE)
  {
    tmp = tmp[-c(1,nrow(tmp)),-c(1,ncol(tmp))]
    tmp1[[1]] = tmp1[[1]][-c(1,length(tmp1[[1]]))]
    tmp2[[1]] = tmp2[[1]][-c(1,length(tmp2[[1]]))]
    
    tmp1[[2]] = tmp1[[2]][-c(1,length(tmp1[[2]]))]
    tmp2[[2]] = tmp2[[2]][-c(1,length(tmp2[[2]]))]
  }
  

 tmp = matrix(as.vector(tmp), ncol = ncol(tmp))
  new("BinnDepth2d", freq = tmp, mid_x =  tmp1[[2]], mid_y = tmp2[[2]], breaks_x = tmp1[[1]], breaks_y = tmp2[[1]], input_data = x,
      max_depth_x = tmp1[[3]], max_depth_y = tmp2[[3]])
  #return(result)
}

#'@docType methods
#'@title 2d Binning plot
#'
#'  @param x object of class BinnDepth2d
#'  @param ... graphical parameters passed to plot
#'  @param alpha alpha value for rgb function
#'  @param bg_col backgroud color
#'  @param add_mid logical. If TRUE centers of binns will be marked.
#'
#'@description Binning 2d
#'@export
#'  @seealso \code{\link{depth}}
#'  
#'  @examples
#'  
#'  tmp = binningDepth2D(x = mvrnorm(100,rep(0,2),diag(2)))
#'  plot(tmp)
#'  @keywords
#'  multivariate
#'  nonparametric
#'  robust
#'  depth function
setMethod("plot", signature = c(x = "BinnDepth2d"), function(x,..., alpha = 0.1, bg_col = "red", add_mid = TRUE){
  
  breaks_y = x@breaks_y
  breaks_x = x@breaks_x
  breaks_x[is.infinite(breaks_x)] = extendrange(x@input_data[,1], f=10)
  breaks_y[is.infinite(breaks_y)] = extendrange(x@input_data[,2], f=10)
  
  
  xlim = extendrange(x@input_data[,1],f=0.1)
  ylim = extendrange(x@input_data[,2],f=0.1)
  plot(x@input_data, xlim = xlim, ylim = ylim, xlab = "", ylab = "")
  rect(xleft=min(breaks_x),xright=max(breaks_x), ybottom=min(breaks_y), ytop = max(breaks_y), col = rgb(1, 0, 0, alpha))
  
  segments(x0=min(breaks_x), x1 = max(breaks_x), y0 = breaks_y, y1= breaks_y, lty = 2)
  segments(x0=breaks_x, x1 = breaks_x, y0 = min(breaks_y), y1= max(breaks_y), lty = 2)
  
  tmp = expand.grid(x@mid_x,x@mid_y)
  if(add_mid) points(tmp[,1], tmp[,2], col = "red", pch = 17)
  #abline(v = breaks_x, lty = 2)
  #abline(h = breaks_y, lty = 2)
})

