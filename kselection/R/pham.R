##
##  Implementation of the 2004 Pham et al. paper
##
##  Created by Daniel Rodriguez Perez on 6/9/2014.
##
##  Copyright (c) 2014 Daniel Rodriguez Perez.
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>
##

#' Selection of K in K-means Clustering
#' 
#' Selection of k in k-means clustering based on Pham et al. paper.
#' 
#' @param x numeric matrix of data, or an object that can be coerced to such a
#'        matrix.
#' @param fun_cluster function to cluster by (e.g. \code{kmeans}). The first 
#'        parameter of the function must a numeric matrix and the second the
#'        number of clusters. The function must return an object with a named
#'        attribute \code{withinss} which is a numeric vector with the within.
#' @param max_centers maximum number of clusters for evaluation.
#' @param k_threshold maximum value of \eqn{f(K)} from which can not be
#'        considered the existence of more than one cluster in the data set.
#'        The default value is 0.85.
#' @param progressBar show a progress bar.
#' @param trace display a trace of the progress.
#' @param parallel If set to true, use parallel \code{foreach} to execute the
#'        function that implements the kmeans algorithm. Must register parallel
#'        before hand, such as \code{doMC} or others. Selecting this option the
#'        progress bar is disabled.
#' @param ... arguments to be passed to the kmeans method.
#' 
#' @return an object with the \eqn{f(K)} results.
#' 
#' @details
#' This function implements the method proposed by Pham, Dimov and Nguyen for
#' selecting the number of clusters for the K-means algorithm. In this method
#' a function \eqn{f(K)} is used to evaluate the quality of the resulting
#' clustering and help decide on the optimal value of \eqn{K} for each data
#' set. The \eqn{f(K)} function is defined as
#' \deqn{f(K) = \left\{
#' \begin{array}{rl}
#'  1 & \mbox{if $K = 1$} \\
#'  \frac{S_K}{\alpha_K S_{K-1}} & \mbox{if $S_{K-1} \ne 0$, $\forall K >1$} \\
#'  1 & \mbox{if $S_{K-1} = 0$, $\forall K >1$}
#' \end{array} \right.}{f(K) =
#'  1, if K = 1;
#'  (S_K)/(\alpha_K S_{K-1}, if S_{K-1} \ne 0, forall K >1;
#'  1, if S_{K-1} = 0, forall K > 1}
#' where \eqn{S_K} is the sum of the distortion of all cluster and \eqn{\alpha_K}
#' is a weight factor which is defined as
#' \deqn{\alpha_K = \left\{
#' \begin{array}{rl}
#'  1 - \frac{3}{4 N_d}                        & \mbox{if $K = 1$ and $N_d > 1$} \\
#'  \alpha_{K-1} + \frac{1 - \alpha_{K-1}}{6}  & \mbox{if $K > 2$ and $N_d > 1$}
#' \end{array} \right.}{\alpha_K = 
#'  1 - 3/(4 * N_d), if K = 1 and N_d > 1;
#'  \alpha_{K-1} + (1 - \alpha_{K-1})/6, if K > 2 and N_d > 1}
#' where \eqn{N_d} is the number of dimensions of the data set.
#' 
#' In this definition \eqn{f(K)} is the ratio of the real distortion to the
#' estimated distortion and decreases when there are areas of concentration in
#' the data distribution.
#' 
#' The values of \eqn{K} that yield \eqn{f(K) < 0.85} can be recommended for
#' clustering. If there is not a value of \eqn{K} which \eqn{f(K) < 0.85}, it
#' cannot be considered the existence of clusters in the data set.
#' 
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'
#' # Execute the method
#' sol <- kselection(dat)
#' 
#' # Get the results
#' k   <- num_clusters(sol) # optimal number of clustes
#' f_k <- get_f_k(sol)      # the f(K) vector
#' 
#' # Plot the results
#' plot(sol)
#' 
#' \dontrun{
#' # Parallel
#' require(doMC)
#' registerDoMC(cores = 4)
#' 
#' system.time(kselection(dat, max_centers = 50 , nstart = 25))
#' system.time(kselection(dat, max_centers = 50 , nstart = 25, parallel = TRUE))
#' }
#' 
#' @author Daniel Rodriguez
#' 
#' @references
#' D T Pham, S S Dimov, and C D Nguyen, "Selection of k in k-means clustering",
#' Mechanical Engineering Science, 2004, pp. 103-119.
#' 
#' @seealso \code{\link{num_clusters}}, \code{\link{get_f_k}}
#' 
#' @import tools
#' 
#' @rdname kselection
#' @export kselection
kselection <- function(x,
                       fun_cluster = stats::kmeans,
                       max_centers = 15,
                       k_threshold = 0.85,
                       progressBar = FALSE,
                       trace       = FALSE,
                       parallel    = FALSE, ...) {
  if (!is.function(fun_cluster))
    stop("'fun_cluster' must be a function.")
  
  if (max_centers < 2)
    stop("'max_centers' must be greater than 2")
  
  if (!is.logical(progressBar)) {
    progressBar <- FALSE
    warning("'progressBar' must be a logical")
  }
  
  if (!is.logical(trace)) {
    trace <- FALSE
    warning("'trace' must be a logical")
  }
  
  if (!is.logical(parallel)) {
    parallel <- FALSE
    warning("'parallel' must be a logical")  
  }
  
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop('x must contain numerical data')
  
  num_row <- dim(x)[1]
  num_col <- dim(x)[2]
  
  if (num_row <= max_centers) {
    warning('The maximum number of clusters has been reduced from ', max_centers, ' to ', num_row - 1)
    max_centers <- num_row - 1
  } 
  
  f_k <- rep(1, max_centers)
  s_k <- rep(1, max_centers)
  a_k <- alpha_k(num_col, max_centers)
    
  if (parallel && requireNamespace('foreach', quietly = TRUE)) {
    s_k <- foreach::"%dopar%"(
      foreach::foreach(k = 1:max_centers, .combine = c), {
        mod_info <- fun_cluster(x, k, ...)
        s_k      <- sum(mod_info$withinss)
      })
  } else {
    if (progressBar) {
      pb <- txtProgressBar(min   = 0,
                           max   = max_centers,
                           style = 3)
    }
    
    for (k in 1:max_centers) {
      mod_info <- fun_cluster(x, k, ...)
      s_k[k]   <- sum(mod_info$withinss)
      
      if (progressBar) {
        setTxtProgressBar(pb, k)
      }
    }
  }
  
  for (k in 1:max_centers) {    
    if (k == 1) {
      f_k[k] <- 1  
    } else {
      if (s_k[k - 1] == 0) {
        f_k[k] <- 1
      } else {
        f_k[k] <- s_k[k]  / (a_k[k] * s_k[k - 1])
      }      
    }
    
    if (trace) {
      message('Number of clusters', k,' - f(k) = ', f_k[k])
    }
  }
  
  result <- list(k           = which_cluster(f_k, k_threshold),
                 f_k         = f_k,
                 max_centers = max_centers,
                 k_threshold = k_threshold,
                 fun_cluster = fun_cluster)
  class(result) <- 'Kselection'
  
  return(result)
}

#' Get the \code{k_threshold}
#'
#' Get the maximum value of \eqn{f(K)} from which can not be considered the
#' existence of more than one cluster.
#'
#' @param obj the output of \code{kselection} function.
#'
#' @return the \code{k_threshold} value.
#' 
#' @author Daniel Rodriguez
#' 
#' @seealso \code{\link{set_k_threshold}}
#'
#' @rdname get_k_threshold
#' @export get_k_threshold
get_k_threshold <- function(obj) {
  UseMethod("get_k_threshold")
}

#' @method get_k_threshold default
#' @export
get_k_threshold.default <- function(obj) {
  NULL
}

#' @method get_k_threshold Kselection
#' @export
get_k_threshold.Kselection <- function(obj) {
  obj$k_threshold
}

#' Set the \code{k_threshold}
#'
#' Set the maximum value of \eqn{f(K)} from which can not be considered the
#' existence of more than one cluster.
#'
#' @param obj the output of \code{kselection} function.
#' @param k_threshold maximum value of \eqn{f(K)} from which can not be
#'        considered the existence of more than one cluster in the data set.
#'
#' @return the output of kselection function with new \code{k_threshold}.
#' 
#' @author Daniel Rodriguez
#' 
#' @seealso \code{\link{get_k_threshold}}
#'
#' @rdname set_k_threshold
#' @export set_k_threshold
set_k_threshold <- function(obj, k_threshold) {
  UseMethod("set_k_threshold")
}

#' @method set_k_threshold Kselection
#' @export
set_k_threshold.Kselection <- function(obj, k_threshold) {
  if (length(k_threshold) != 1L)
    stop('k_threshold must be scalar')
  if (!is.numeric(k_threshold))
    stop('k_threshold must be numeric')
  if (k_threshold <= 0)
    stop('k_threshold must be numeric bigger than 0')

  obj$k_threshold <- k_threshold
  obj
}

#' Get the \eqn{f(K)} vector
#'
#' Get the \eqn{f(K)} vector.
#'
#' @param obj the output of \code{kselection} function.
#'
#' @return the vector of \eqn{f(K)} function.
#'
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'
#' # Get the f(k) vector
#' sol <- kselection(dat)
#' f_k <- get_f_k(sol)
#' 
#' @author Daniel Rodriguez
#' 
#' @seealso \code{\link{num_clusters}}, \code{\link{num_clusters_all}}
#'
#' @rdname get_f_k
#' @export get_f_k
get_f_k <- function(obj) {
  UseMethod("get_f_k")
}

#' @method get_f_k default
#' @export
get_f_k.default <- function(obj) {
  NULL
}

#' @method get_f_k Kselection
#' @export
get_f_k.Kselection <- function(obj) {
  obj$f_k
}

#' Get the number of clusters.
#'
#' The optimal number of clusters proposed by the method.
#'
#' @param obj the output of kselection function.
#'
#' @return the number of clusters proposed.
#'
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'
#' # Get the optimal number of clustes
#' sol <- kselection(dat)
#' k   <- num_clusters(sol)
#' 
#' @author Daniel Rodriguez
#' 
#' @seealso \code{\link{num_clusters_all}}, \code{\link{get_f_k}}
#'
#' @rdname num_clusters
#' @export num_clusters
num_clusters <- function(obj) {
  UseMethod("num_clusters")
}

#' @method num_clusters default
#' @export
num_clusters.default <- function(obj) {
  NULL
}

#' @method num_clusters Kselection
#' @export
num_clusters.Kselection <- function(obj) {
  obj$k
}

#' Get all recommended numbers of clusters
#'
#' The number of cluster which could be recommender according the method
#' threshold.
#'
#' @param obj the output of \code{kselection} function.
#'
#' @return an array of number of clusters that could be recommended.
#'
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'
#' # Get the optimal number of clustes
#' sol <- kselection(dat)
#' k   <- num_clusters(sol)
#' 
#' @author Daniel Rodriguez
#' 
#' @seealso \code{\link{num_clusters}}, \code{\link{get_f_k}}
#' 
#' @rdname num_clusters_all
#' @export num_clusters_all
num_clusters_all <- function(obj) {
  UseMethod("num_clusters_all")
}

#' @method num_clusters_all default
#' @export
num_clusters_all.default <- function(obj) {
  NULL
}

#' @method num_clusters_all Kselection
#' @export
num_clusters_all.Kselection <- function(obj) {
  which(get_f_k(obj) < get_k_threshold(obj))
}

#' @method plot Kselection
#' @export
plot.Kselection <- function(x, ...) {
  max_y   <- 1.1 * max(x$f_k)
  valid_k <- num_clusters_all(x)

  plot(x$f_k,
       main = to.string(x),
       type = 'b',
       xlab = 'Number of clusters k',
       ylab = 'f(K)',
       ylim = c(0, max_y))
  lines(x$k, x$f_k[x$k],
        col  = 'green',
        pch  = 19,
        type = 'p')
  if (length(valid_k) > 0) {
    lines(valid_k, x$f_k[valid_k],
          col  = 'green',
          type = 'p')
    if (x$f_k[length(x$f_k)] > max_y / 2)
      legend('bottomright', c('Lower f(K)', 'Recommended K'),
             col = c('green','green'),
             pch = c(19, 1))
    else
      legend('topright', c('Lower f(K)', 'Recommended K'),
             col = c('green','green'),
             pch = c(19, 1))
  }
}

#' @method print Kselection
#' @export
print.Kselection <- function(x, ...) {
  cat(to.string(x))
}

# Calculate the weight factor for kselection
#
# Calculate the weight factor for kselection.
#
# @param n_d the number of dimensions (attributes).
# @param k the number of iterations.
#
# @return an array with the weights.
# 
# @author Daniel Rodriguez
alpha_k <- function(n_d, k) {
  result <- 1 - 3/(4 * n_d)
  result <- rep(result, k)
  
  for (i in 3:k) {
    result[i] <- result[i - 1] + (1 - result[i - 1]) / 6
  }

  return(result)
}

# Calculate the optimal cluster for kselection
#
# Calculate the optimal cluster for kselection.
#
# @param f_k the \eqn{f(K)} array.
# @param k_threshold cluster selection threshold for \eqn{f(K)} value.
#
# @return the number of clusters.
# 
# @author Daniel Rodriguez
which_cluster <- function(f_k, k_threshold) {
  k <- which(f_k == min(f_k) & f_k < k_threshold)
  if (length(k) == 0)
    return(1)
  else
    return(min(k))
}

# Generate a string with the recomender number of k
#
# Generate a string with the recomender number of k.
#
# @param obj the output of kselection function.
#
# @return an string with the recomendation.
#
# @author Daniel Rodriguez
to.string <- function(obj) {
  if (num_clusters(obj) == 1)
    paste('f(k) finds', num_clusters(obj), 'cluster')
  else
    paste('f(k) finds', num_clusters(obj), 'clusters')
}
