# #### 
# #### Computes the Takens' estimator.
# ####   
# #### details
# #### The Takens' estimator is an alternative to determining a finite correlation dimension suggested by Takens.
# #### Takens' estimator provides a maximum likelihood estimate and it is also used as statistic for nonlinearity
# #### testing.
# #### references Kantz, H.  and Schreiber, T.: Nonlinear Time series Analysis (Cambridge university press)
# #### Theiler, J. Lacunarity in a best estimator of fractal dimension. Phys. Lett. A,135, 195 (1988)
# ####
# #### author Constantino A. Garcia
# 
# takens.estimator = function(time.series, min.embedding.dim = 2, max.embedding.dim = 5, time.lag = 1,
#                             min.radius, max.radius, n.points.radius, do.plot = TRUE, 
#                             theiler.window = 0, number.boxes = NULL){
#   
#   
#     corr.dim = corrDim(time.series=time.series,min.embedding.dim= min.embedding.dim,
#                        max.embedding.dim=max.embedding.dim, time.lag = time.lag,
#                         min.radius = min.radius, max.radius = max.radius,
#                         n.points.radius = n.points.radius, do.plot=do.plot,
#                         theiler.window = theiler.dist, number.boxes =number.boxes)
# 
#     # vector that will store the takens' estimator for each dimension
#     takens.estimator = vector(mode = "numeric", length = nrow(corr.dim) )
#     names(takens.estimator) = min.embedding.dim:max.embedding.dim
#     # the radius axis. We have to revert it since it's sorted in decreasing order
#     radius.vector =  rev(as.numeric(dimnames(corr.dim)[[2]]))
#     # for each dimension of the correlation dimension matrix, we will calculate
#     # the takens' estimator
#     position = 1
#     for (dimension in min.embedding.dim:max.embedding.dim){
#       current.dimension = as.character(dimension)
#       integrand = rev(corr.dim[current.dimension,])/radius.vector
#       integral = trapezoidalRule(c(0,radius.vector), c(0,integrand))
#       lc=log(corr.dim[current.dimension,])
#       le=log(rev(radius.vector))
#       cat("*\n")
#       takens.estimator[[current.dimension]] = corr.dim[current.dimension,1]/integral
#     }
#     print(takens.estimator)
#     return (list(takens.estimator=takens.estimator,correlation.matrix = corr.dim))
# }



#' Time Reversibility statistic
#' @details  The time simmetry statistic  measures the asymmetry of a time series 
#' under time reversal by implementing the third order statistic:
#' \deqn{E[s_n - s_{n-\tau})^3] }{E [s(n) - s(n-tau))^3]}.
#' Since linear stochastic series are symmetric under time reversal, this statistic
#' may be used for testing the assertion that the data was generated
#' from a stationary linear stochastic process or not.
#' @param time.series The time series used to compute the statistic
#' @param tau Time delay used to compute the third order statistic.
#' @return The time simmetry statistic for the delays specified with
#' \emph{tau}.
#' @seealso \code{\link{timeAsymmetry}}
#' @author Constantino A. Garcia
#' @rdname timeAsymmetry2
#' @export timeAsymmetry2
timeAsymmetry2 <- function(time.series,tau){
  len.ts = length(time.series)
  statistic = c()
  for (i in seq_along(tau)){
    statistic[[i]] = mean( ( time.series[(tau[[i]]+1):len.ts] - time.series[1:(len.ts-tau[[i]])] )^3)  
  }
  return(statistic)
}

#' Time Reversibility statistic
#' @details  The time simmetry statistic  measures the asymmetry of a time series 
#' under time reversal by calculating:
#' \deqn{E[s_n\cdot s_{n+1}^2] - E[s_n^2\cdot s_{n+1}] }{E[s_n * s_{n+1}^2] - E[s_n^2 * s_{n+1}] }.
#' Since linear stochastic series are symmetric under time reversal, this statistic
#' may be used for testing the assertion that the data was generated
#' from a stationary linear stochastic process or not.
#' @param time.series The time series used to compute the statistic.
#' @return The time simmetry statistic.
#' @author Constantino A. Garcia
#' @rdname timeAsymmetry
#' @export timeAsymmetry
timeAsymmetry = function(time.series){
  len=length(time.series)
  mean(time.series[1:(len-1)] *
         time.series[2:len]^2 -  time.series[1:(len-1)]^2 * time.series[2:len])
}