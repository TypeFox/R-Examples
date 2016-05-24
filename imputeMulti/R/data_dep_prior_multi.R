#' @title Data Dependent Prior for Multinomial Distribution
#' @description Creates a data depedent prior for p-dimensional multinomial distributions
#' using a conjugate prior (eg \eqn{Dirichlet(\alpha)}) based on 20% of the data.
#' @param dat A \code{data.frame}. All variables must be factors
#' @references Darnieder, William Francis. Bayesian methods for data-dependent priors. 
#' Dissertation. The Ohio State University, 2011. 
#' @return A \code{data.frame} containing identifiers for all possible \eqn{P(Y=y)} and 
#' the associated prior-counts, \eqn{\alpha}
#' @seealso \code{\link{expand.grid}}
#' @export
data_dep_prior_multi <- function(dat) {
  if (!all(apply(dat, 2, is.factor))) {
    # enforce factor variables
    dat <- data.frame(apply(dat, 2, function(x) as.factor(x)))
  }
  
  enum <- expand.grid(sapply(dat, levels))
  
  comp <- which(stats::complete.cases(dat))
  comp_frac <- length(comp) / nrow(dat) 
  
  if (comp_frac < .2) {
    prior <- count_levels(dat, enum_list= enum, hasNA= "no") 
  } else {
    n <- round(.2 * length(comp) / comp_frac, 0)
    samp <- sample(comp, size= n)
    prior <- count_levels(dat[samp,], enum_list= enum, hasNA= "no")
  }
  
  prior <- merge(enum, prior, all.x= TRUE, all.y=FALSE)
  names(prior)[ncol(prior)] <- "alpha" # naming convention of dirichlet prior
  prior$alpha <- ifelse(is.na(prior$alpha), 1, prior$alpha)
  
  return(prior)
}


# helper function for checking priors in 3 main functions. 
# reduces duplication of code
check_prior <- function(dat, conj_prior= c("none", "data.dep", "flat.prior", "non.informative"),
                        alpha= NULL, verbose= FALSE,
                        outer= FALSE, enum_comp= NULL) {
  if (outer) { # called w/in multinomial_impute
    if (conj_prior == "none") return(NULL)
    
    if (conj_prior == "data.dep") {
      if (!is.null(alpha)) {
        if (verbose == TRUE) print("Using user-supplied data dependent prior.")
        message("Using user-supplied data dependent prior.")
      } else {
        if (verbose == TRUE) print("Calculating data dependent prior.")
        message("Calculating data dependent prior.")
        alpha <- data_dep_prior_multi(dat= dat)
      }
    } else if (conj_prior == "flat.prior") {
      if (!(is.vector(alpha) & length(alpha) == 1)) {
        stop("Flat priors must be supplied as a scalar.")
      }
      alpha <- alpha
    } else if (conj_prior == "non.informative") {
      alpha <- 1 # Jeffrey's prior * 2
    }
    return(alpha)
  } else { # called w/in multinomial_em or multinomial_data_aug
    if (conj_prior != "none") {
      if (conj_prior == "data.dep") {
        if (nrow(alpha) != nrow(enum_comp)) {
          stop("nrow(alpha) must match nrow(enum_comp).")
        }
        enum_comp <- merge(enum_comp, alpha)
      } else if (conj_prior == "flat.prior") {
        if (!(is.vector(alpha) & length(alpha) == 1)) {
          stop("Flat priors must be supplied as a scalar.")
        }
        enum_comp$alpha <- alpha
      } else if (conj_prior == "non.informative") {
        enum_comp$alpha <- 1 # Jeffrey's prior * 2
      } 
      # calc theta_y from alpha
      enum_comp$theta_y <- enum_comp$alpha / sum(enum_comp$alpha)
    } else {
      enum_comp$theta_y <- stats::runif(nrow(enum_comp))
      enum_comp$theta_y <- enum_comp$theta_y / sum(enum_comp$theta_y)
    }
    return(enum_comp)
  }
}