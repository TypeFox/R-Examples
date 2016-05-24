###################################################################### 
#' @title Bootstrapped multi-way standard error clustering
#' 
#' @description Return a bootstrapped multi-way cluster-robust variance-covariance matrix
#'
#' @param model The estimated model, usually an \code{lm} or \code{glm} class object
#' @param cluster A vector, \code{matrix}, or \code{data.frame} of cluster variables,
#' where each column is a separate variable.  If the vector \code{1:nrow(data)}
#' is used, the function effectively produces a regular 
#' heteroskedasticity-robust matrix.
#' @param parallel Scalar or list.  If a list, use the list as a list
#' of connected processing cores/clusters.  Scalar values of \code{TRUE}
#' and \code{"snow"} (which are equivalent) ask \code{boot} to handle parallelization, as does
#' \code{"multicore"}.  See the \code{parallel} and \code{boot} package.
#' @param use_white Logical or \code{NULL}.  See description below.
#' @param force_posdef Logical.  Force the eigenvalues of the variance-covariance
#' matrix to be positive.
#' @param R \code{Integer}. The number of bootstrap replicates; passed directly to \code{boot}.
#' @param boot_type \code{"xy"}, \code{"residual"}, or \code{"wild"}.  See details.
#' @param wild_type \code{"rademacher"}, \code{"mammen"}, or \code{"norm"}.  See details.
#' @param debug Logical.  Print internal values useful for debugging to 
#' the console.
#'
#' @keywords clustering multi-way robust standard errors bootstrap boot block
#' 
#' @details
#' This function implements cluster bootstrapping (also known as the block bootstrap)
#' for variance-covariance matrices, following Cameron, Gelbach, & Miller (CGM) (2008).
#' Usage is generally similar to the \code{cluster.vcov} function in this package, but this
#' function does not support degrees of freedome corrections or leverage adjustments.
#' 
#' In the terminology that CGM (2008) use, this function implements
#' \emph{pairs, residual, or wild cluster bootstrap-se}.
#' 
#' A pairs (or xy) cluster bootstrap can be obtained by setting \code{boot_type = "xy"}, 
#' which resamples the entire regression data set (both X and y).
#' Setting \code{boot_type = "residual"} will obtain a residual cluster
#' bootstrap, which resamples only the residuals (in this case, we resample the blocks/clusters
#' rather than the individual observations' residuals).  To get a wild cluster bootstrap set
#' \code{boot_type = "wild"}, which does not resample anything, but instead reforms the
#' dependent variable by multiplying the residual by a randomly drawn value and adding the
#' result to the fitted value.  The default method is the pairs/xy bootstrap.
#' 
#' There are three built-in distributions to draw multipliers from for wild bootstraps:
#' the Rademacher (\code{wild_type = "rademacher"}, the default), which draws from [-1, 1],
#' each with P = 0.5, Mammen's suggested distribution (\code{wild_type = "mammen"}, see 
#' Mammen, 1993), and the standard normal/Gaussian distribution (\code{wild_type = "norm"}).
#' The default is the Rademacher distribution, following CGM (2008).  Alternatively, you can
#' set the function to draw multipliers from by assigning
#' \code{wild_type} to a function that takes no arguments and returns a single real value.
#' 
#' Multi-way clustering is handled as described by Petersen (2009) and generalized
#' according to Cameron, Gelbach, & Miller (2011).  This means that cluster.boot
#' estimates a set of variance-covariance matrices \emph{for the variables} separately 
#' and then sums them (subtracting some matrices and adding others).
#' The method described by CGM (2011) estimates a set of variance-covariance matrices
#' \emph{for the residuals} (sometimes referred to as the meat of the sandwich estimator)
#' and sums them appropriately.  Whether you sum the meat matrices and then compute
#' the model's variance-covariance matrix or you compute a series of model matrices
#' and sum those is mathematically irrelevant, but may lead to (very) minor numerical
#' differences.
#' 
#' Instead of passing in a vector, matrix, data.frame, etc, to specify the cluster variables,
#' you can use a formula to specify which variables from the
#' original data frame to use as cluster variables, e.g., \code{~ firmid + year}.
#' 
#' Ma (2014) suggests using the White (1980) 
#' variance-covariance matrix as the final, subtracted matrix when the union 
#' of the clustering dimensions U results in a single observation per group in U;
#' e.g., if clustering by firm and year, there is only one observation
#' per firm-year, we subtract the White (1980) HC0 variance-covariance
#' from the sum of the firm and year vcov matrices.  This is detected
#' automatically (if \code{use_white = NULL}), but you can force this one way 
#' or the other by setting \code{use_white = TRUE} or \code{FALSE}.
#' 
#' Unlike the \code{cluster.vcov} function, this function does not depend upon the 
#' \code{\link[sandwich]{estfun}}
#' function from the \CRANpkg{sandwich} package, although it does make use of the \code{\link[sandwich]{vcovHC}}
#' function for computing White (1980) variance-covariance matrices. 
#' 
#' Parallelization (if used) is handled by the \bold{boot} package.  Be sure to set
#' \code{options(boot.ncpus = N)} where \code{N} is the number of CPU cores you want
#' the \code{boot} function to use.
#' 
#' @return a \eqn{K x K} variance-covariance matrix of type \code{matrix}
#' 
#' @export
#' 
#' @seealso \code{\link{cluster.vcov}} for clustering using asymptotics 
#' 
#' @author Nathaniel Graham \email{npgraham1@@gmail.com}
#' 
#' @references
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based improvements 
#' for inference with clustered errors. The Review of Economics and Statistics, 90(3), 414-427.
#' \doi{10.1162/rest.90.3.414}
#' 
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). Robust inference with multiway 
#' clustering. Journal of Business & Economic Statistics, 29(2).
#' \doi{10.1198/jbes.2010.07136}
#' 
#' Mammen, E. (1993). Bootstrap and wild bootstrap for high dimensional linear models. The 
#' Annals of Statistics, 255-285.
#' \doi{10.1214/aos/1176349025}
#' 
#' Petersen, M. A. (2009). Estimating standard errors in finance panel data sets: Comparing 
#' approaches. Review of Financial Studies, 22(1), 435-480.
#' \doi{10.1093/rfs/hhn053}
#' 
#' White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct 
#' test for heteroskedasticity. Econometrica: Journal of the Econometric Society, 817--838.
#' \doi{10.2307/1912934}
#' 
#' @importFrom sandwich vcovHC vcovHC.default
#' @importFrom boot boot
#' @importFrom compiler cmpfun
#' @importFrom parallel clusterExport
#' @importFrom stats coef cov model.frame residuals formula rnorm fitted update.formula na.pass
#' @importFrom utils combn
#' 
#' @examples
#' \dontrun{
#' library(lmtest)
#' data(petersen)
#' m1 <- lm(y ~ x, data = petersen)
#' 
#' # Cluster by firm
#' boot_firm <- cluster.boot(m1, petersen$firmid)
#' coeftest(m1, boot_firm)
#' 
#' # Cluster by firm using a formula
#' boot_firm <- cluster.boot(m1, ~ firmid)
#' coeftest(m1, boot_firm)
#' 
#' # Cluster by year
#' boot_year <- cluster.boot(m1, petersen$year)
#' coeftest(m1, boot_year)
#' 
#' # Double cluster by firm and year
#' boot_both <- cluster.boot(m1, cbind(petersen$firmid, petersen$year))
#' coeftest(m1, boot_both)
#' 
#' # Cluster by firm with wild bootstrap and custom wild distribution
#' boot_firm2 <- cluster.boot(m1, petersen$firmid, boot_type = "wild",
#'                            wild_type = function() sample(c(-1, 1), 1))
#' coeftest(m1, boot_firm)
#' 
#' # Go multicore using the parallel package
#' require(parallel)
#' cl <- makeCluster(4)
#' options(boot.ncpus = 4)
#' boot_both <- cluster.boot(m1, cbind(petersen$firmid, petersen$year), parallel = cl)
#' stopCluster(cl)
#' coeftest(m1, boot_both)
#' 
#' # Go multicore using the parallel package, let boot handle the parallelization
#' require(parallel)
#' options(boot.ncpus = 8)
#' boot_both <- cluster.boot(m1, cbind(petersen$firmid, petersen$year), parallel = TRUE)
#' coeftest(m1, boot_both)
#' }
cluster.boot <- function(model, cluster, parallel = FALSE, use_white = NULL, 
                         force_posdef = FALSE, R = 300, boot_type = "xy", wild_type = "rademacher",
                         debug = FALSE) {
  
  if(inherits(cluster, "formula")) {
    cluster_tmp <- expand.model.frame(model, cluster, na.expand = FALSE)
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  } else {
    cluster <- as.data.frame(cluster, stringsAsFactors = FALSE)
  }
  
  cluster_dims <- ncol(cluster)
  
  # total cluster combinations, 2^D - 1
  tcc <- 2 ** cluster_dims - 1
  # all cluster combinations
  acc <- list()
  for(i in 1:cluster_dims) {
    acc <- append(acc, combn(1:cluster_dims, i, simplify = FALSE))
  }
  if(debug) print(acc)
  
  # We need to subtract matrices with an even number of combinations and add
  # matrices with an odd number of combinations
  vcov_sign <- sapply(acc, function(i) (-1) ** (length(i) + 1))
  
  # Drop the original cluster vars from the combinations list
  acc <- acc[-1:-cluster_dims]
  if(debug) print(acc)
  
  # Handle omitted or excluded observations
  if(!is.null(model$na.action)) {
    if(class(model$na.action) == "exclude") {
      cluster <- cluster[-model$na.action,]
    } else if(class(model$na.action) == "omit") {
      cluster <- cluster[-model$na.action,]
    }
    cluster <- as.data.frame(cluster)  # silly error somewhere
  }
  if(debug) print(class(cluster))
  
  # Factors in our clustering variables can potentially cause problems
  # Blunt fix is to force conversion to characters
  i <- !sapply(cluster, is.numeric)
  cluster[i] <- lapply(cluster[i], as.character)
  
  # Make all combinations of cluster dimensions
  if(cluster_dims > 1) {
    for(i in acc) {
      cluster <- cbind(cluster, Reduce(paste0, cluster[,i]))
    }
  }
  
  df <- data.frame(M = integer(tcc),
                   N = integer(tcc),
                   K = integer(tcc))
  
  for(i in 1:tcc) {
    df[i, "M"] <- length(unique(cluster[,i]))
    df[i, "N"] <- length(cluster[,i])
    df[i, "K"] <- model$rank
  }
  
  if(is.null(use_white)) {
    if(cluster_dims > 1 && df[tcc, "M"] == prod(df[-tcc, "M"])) {
      use_white <- TRUE
    } else {
      use_white <- FALSE
    }
  }
  
  if(use_white) {
    df <- df[-tcc,]
    tcc <- tcc - 1
  }
  
  if(debug) {
    print(acc)
    print(paste("Original Cluster Dimensions", cluster_dims))
    print(paste("Theoretical Cluster Combinations", tcc))
    print(paste("Use White VCOV for final matrix?", use_white))
  }
  

  if("model" %in% names(model)) {
    full_data <- model$model
  } else {
    full_data <- model.frame(model)
  }
  
  boot.outs <- list()
  
  par_type <- "no"
  if(length(parallel) > 1) {
    par_cluster <- parallel
    clusterExport(par_cluster, varlist = c("cluster", "model"), envir = environment())
    par_type <- "snow"
  } else if(parallel == TRUE || parallel == "snow") {
    par_type <- "snow"
    par_cluster <- NULL
  } else if(parallel == "multicore") {
    par_type <- "multicore"
    par_cluster <- NULL
  }
  
  dat_loc <- which(names(model$call) == "data") 
  args <- model$call[c(-1, -dat_loc)]
  wild_func <- function() NULL
  
  # Setup the bootstrap functions, depending upon the type of bootstrap requested
  if(boot_type == "xy") {
    est.func <- cmpfun(function(grp, i, data2, clustvar, reg_arglist, boot_args) {
      j <- unlist(lapply(i, function(n) which(n == clustvar)))
      coef(boot_args$estimator(reg_arglist, data = data2[j,]))
    })
  } else if(boot_type == "residual") {
    args$formula <- update.formula(formula(model), y_boot ~ .)
    full_data$y_boot <- 0
    
    est.func <- cmpfun(function(grp, i, data2, clustvar, reg_arglist, boot_args) {
      j <- unlist(lapply(i, function(n) which(n == clustvar)))
      data2$y_boot <- fitted(boot_args$model) + residuals(boot_args$model)[j]
      coef(boot_args$estimator(reg_arglist, data = data2))
    })
  } else if(boot_type == "wild") {
    args$formula <- update.formula(formula(model), y_boot ~ .)
    full_data$y_boot <- 0
    
    if(is.function(wild_type)) {
      wild_func <- wild_type
    } else {
      if(wild_type == "rademacher") {
        wild_func <- cmpfun(function() sample(c(-1, 1), 1))
      } else if(wild_type == "mammen") {
        wild_func <- cmpfun(function() sample(c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2), 1,
                                              prob = c((sqrt(5) + 1) / (2 * sqrt(5)), 
                                                       (sqrt(5) - 1) / (2 * sqrt(5)))))
      } else if(wild_type == "norm") {
        wild_func <- cmpfun(function() rnorm(1))
      }
    }
    
    
    est.func <- cmpfun(function(grp, i, data2, clustvar, reg_arglist, boot_args) {
      j <- unlist(lapply(grp, function(n) rep_len(boot_args$wild_func(), sum(n == clustvar))))
      data2$y_boot <- fitted(model) + residuals(model) * j
      coef(boot_args$estimator(reg_arglist, data = data2))
    })
  }
  boot_args <- new.env()
  boot_args$estimator <- eval(model$call[[1]])
  boot_args$model <- model
  boot_args$wild_func <- wild_func
  
  for(i in 1:tcc) {
    boot.outs[[i]] <- boot(unique(cluster[,i]), est.func, R = R,
                           parallel = par_type, cl = par_cluster,
                           data2 = full_data, clustvar = cluster[,i],
                           reg_arglist = args, boot_args = boot_args)
  }
  
  if(debug) {
    print(df)
    print(vcov_sign)
  }
  
  vcov_matrices <- list()
  for(i in 1:tcc) {
    vcov_matrices[[i]] <- vcov_sign[i] * cov(boot.outs[[i]]$t)
  }
  
  if(use_white) {
    i <- i + 1
    vcov_matrices[[i]] <- vcov_sign[i] * vcovHC(model, type = "HC0")
  }
  
  if(debug) {
    print(vcov_matrices)
  }
  
  vcov_matrix <- Reduce('+', vcov_matrices)
  
  if(force_posdef) {
    decomp <- eigen(vcov_matrix, symmetric = TRUE)
    if(debug) print(decomp$values)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    vcov_matrix <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }
  
  return(vcov_matrix)
}