#' Generator functions
#' 
#' These functions are intended to be used with \code{\link{sim_gen}} and not
#' interactively. They are designed to draw random numbers according to the
#' setting of grouping variables.
#'  
#' @param mean the mean passed to the random number generator, for example
#'   \code{\link{rnorm}}.
#' @param sd the standard deviation passed to the random number generator, for
#'   example \link{rnorm}.
#' @param name name of variable as character in which random numbers are stored.
#' @param rho the correlation used to create the variance covariance matrix for
#'   a SAR process - see \code{\link[spdep]{cell2nb}}.
#' @param type either "rook" or "queen". See \code{\link[spdep]{cell2nb}} for
#'   details.
#' @param generator a function producing random numbers.
#' @param ... arguments passed to \code{generator}.
#' @param groupVars names of variables as character. Identify groups within
#'   random numbers are constant.
#' @param groupVar a variable name identifying groups.
#' @param timeVar a variable name identifying repeated measurements.
#' 
#' @details \code{gen_norm} is used to draw random numbers from a normal
#'   distribution where all generated numbers are independent.
#' 
#' \code{gen_v_norm} and \code{gen_v_sar} will create an area-level random
#' component. In the case of \code{v_norm}, the error component will be from a
#' normal distribution and i.i.d. from an area-level perspective (all units in
#' an area will have the same value, all areas are independent). v_sar will also
#' be from a normal distribution, but the errors are correlated. The variance
#' covariance matrix is constructed for a SAR(1) - spatial/simultanous
#' autoregressive process. \link[MASS]{mvrnorm} is used for the random number
#' generation. \code{gen_v_norm} and \code{gen_v_sar} expect a variable
#' \code{idD} in the data identifying the areas.
#' 
#' \code{gen_generic} can be used if your world is not normal. You can specify
#' 'any' function as generator, like \code{\link{rnorm}}. Arguments in
#' \code{...} are matched by name or position. The first argument of
#' \code{generator} is expected to be the number of random numbers (not
#' necessarily named \code{n}) and need not to be specified.
#'  
#' @seealso \code{\link{sim_gen}}, \code{\link{sim_gen_x}},
#'   \code{\link{sim_gen_e}}, \code{\link{sim_gen_ec}}, \code{\link{sim_gen_v}},
#'   \code{\link{sim_gen_vc}}, \code{\link[spdep]{cell2nb}}
#'  
#' @rdname generators
#' @export
#'  
#' @examples
#' sim_base() %>% sim_gen_x() %>% sim_gen_e() %>% sim_gen_v() %>% sim_gen(gen_v_sar(name = "vSP"))
#' 
#' # Generic interface
#' set.seed(1)
#' dat1 <- sim(base_id() %>%
#'   sim_gen(gen_generic(rnorm, mean = 0, sd = 4, name = "e")))
#' set.seed(1)
#' dat2 <- sim(base_id() %>% sim_gen_e())
#' all.equal(dat1, dat2)
gen_norm <- function(mean = 0, sd = 1, name = "e") {
  force(mean); force(sd); force(name)
  function(dat) {
    dat <- add_var(dat, rnorm(nrow(dat), mean = mean, sd = sd), name)
    dat
  }
}

add_var <- function(dat, value, name) {
  dat[name] <- value + if(exists(name, dat)) dat[[name]] else 0
  dat
}

#' @rdname generators
#' @export
gen_v_norm <- function(mean = 0, sd = 1, name = "v") {
  force(mean); force(sd); force(name)
  function(dat) {
    tmp <- rnorm(length(unique(dat$idD)), mean = mean, sd = sd)
    dat <- add_var(dat, tmp[dat$idD], name)
    dat
  }
}

#' @rdname generators
#' @export
gen_v_sar <- function(mean = 0, sd = 1, rho = 0.5, type = "rook", name) {
  force(mean); force(sd); force(name); force(rho); force(type)
  function(dat) {
    nDomains <- length(unique(dat$idD))
    # Spatial Structure:
    W <- nb2mat(cell2nb(nDomains, 1, type), style = "W")
    identity <- diag(1, nDomains, nDomains)
    tmp <- identity - rho * W
    sp_var <- sd^2 * chol2inv(chol(crossprod(tmp, tmp)))
    
    # Drawing the numbers:
    v <- mvrnorm(1, mu = rep(mean, length.out = nDomains), Sigma = sp_var)
    
    dat <- add_var(dat, v[dat$idD], name)
    
    dat
  }
}

#' @rdname generators
#' @export
gen_v_ar1 <- function(mean = 0, sd = 1, rho = 0.5, 
                     groupVar = "idD", timeVar = "idT", name) {
  force(mean); force(sd); force(name); force(rho); force(groupVar)
  function(dat) {
    nDomains <- nrow(unique(dat[groupVar]))
    nTime <- length(unique(dat[[timeVar]]))

    # Temporal Structure:
    ar_var <- matrix(0, nTime, nTime)
    for (i in 1:nTime) {
      ar_var[i, ] <- c(rep(0, length.out = (i - 1)), rho^(0:(nTime - i)))
    }
    ar_var <- ar_var + t(ar_var)
    diag(ar_var) <- 1
    ar_var <- 1 / (1 - rho^2) * sd^2 * ar_var
    
    # Drawing the numbers - I assume a completely balanced design:
    dat <- split(dat, dat[groupVar], drop = TRUE) %>%
      lapply(function(df) {
        add_var(
          df, 
          mvrnorm(1, mu = rep(mean, length.out = nTime), Sigma = ar_var), 
          name
        )
      }) %>%
      do.call(what = rbind)
    rownames(dat) <- NULL
    dat[do.call(order, dat[groupVar]), ]
  }
}

#' @rdname generators
#' @export
gen_generic <- function(generator, ..., groupVars = NULL, name) {
  genArgs <- list(...)
  force(generator)
  force(groupVars)
  force(name)
      
  gen_no_group <- function(dat) {
    randomNumbers <- do.call(generator, c(nrow(dat), genArgs))
    add_var(dat, randomNumbers, name)
  }
  
  gen_constant_within_group <- function(dat) {
    dat <- dat %>% arrange_(groupVars)
    nrowGroup <- group_by_(dat, groupVars) %>% group_size
    randomNumbers <- do.call(generator, c(length(nrowGroup), genArgs))
    add_var(dat, rep(randomNumbers, times = nrowGroup), name)
  }
  
  if(is.null(groupVars)) {
    gen_no_group
  } else {
    gen_constant_within_group
  }
}


