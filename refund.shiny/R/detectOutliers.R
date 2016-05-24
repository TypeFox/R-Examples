


#' internal function from 'fda' package
#' 
#' function used in method for fast modified band depth (MBD) calculation
#' 
#' @param n number of columns in your dataset
#' @param p number of rows in your dataset
#' 
#' @author Ying Sun and Marc G.Genton
#'

combinat=function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}


#' fast modified band depth calculation for fda
#'  
#' Method for fast modified band depth (fMBD) calculation
#' 
#' @param data name of dataset
#' 
#' @author Ying Sun and Marc G.Genton
#'

fMBD=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}


#' Identifies outliers for plot_shiny.fosr()
#' 
#' Internal method that assigns band depth values to curves based on exact fast MBD computation (Sun & Genton, 2012).
#' Code modified from fbplot in fda package. 
#' A dataframe of residuals is passed as an argument, and depths and outlying curves are returned
#' 
#' @param data matrix or df of functional observations
#' @param factor a constant that determines the fences for outliers. Defaults to 1.5, as in classical definition for Tukey outliers. 
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @references Sun, Ying, Marc G. Genton, and Douglas W. Nychka. (2012).
#' Exact fast computation of band depth for large functional datasets: How quickly can one million curves be ranked? \emph{Stat}, 1, 68-74.
#' 
#' Sun, Ying, and Marc G. Genton. (2011). Functional boxplots. \emph{Journal of Computational and Graphical Statistics}, 20, 313-334.
#' 

outliers = function(data, factor=1.5){
  # tranpose data so that each column is a curve rather than each row
  data = t(data)
  
  depth = fMBD(data)
  
  tp = dim(data)[1] # number of observations per curve
  n = dim(data)[2]  # number of curves
  x = 1:tp
  
  dp_s = sort(depth, decreasing = TRUE)
  index = order(depth, decreasing = TRUE)
  
  ##### median curve
  med = depth == max(depth)
  medavg = matrix(data[, med], ncol = sum(med), nrow = tp) # gets median curve
  y = apply(medavg, 1, mean)
  ##
  
  ###### get 50% region (analogous to IQR)
  m = ceiling(n * 0.5)
  center = data[, index[1:m]] # 50% region (deepest curves), 'IQR' 
  out = data[, index[(m + 1):n]] # curves outside of 'IQR'
  inf = apply(center, 1, min) 
  sup = apply(center, 1, max)
  ##
  
  ##### get Outliers
  dist = factor * (sup - inf) # defines what it means to be an outlier
  upper = sup + dist
  lower = inf - dist
  outly = (data <= lower) + (data >= upper) # sets matrix, each point in each curve is checked for outlying-ness
  outcol = colSums(outly)
  remove = (outcol > 0)
  colum = 1:n  # lists number of curves 
  outpoint = colum[remove == 1] # gets index of outlying curves
  outcurves = data[, remove] # gets values for outlying curve
  medcurve = data[,med]
  ##
  
 
  # returns index of outliers, subset of data with just outlying curves
  return(list(depth = depth, outpoint = outpoint, medcurve = t(medcurve), outcurves=t(outcurves)))
  
}


