#' @title Generating OpenBUGS Model File
#' 
#' @aliases bugsModel
#' @description Generating OpenBUGS model file
#' 
#' @param formula A formula object, with the DV on the left of an ~ operator, and predictors on the right. For the part on the right of '~', the specification of the submodels can be separated by '|'. So \code{y ~ X1 | X2} means the DV is \code{y},\code{X1} is the term in the mean submodel, and \code{X2} is the term in the dispersion submodel. 
#' @param fd A string that specifies the parent distribution (see \link{cdfqrFamily}).
#' @param sd A string that specifies the sub-family distribution.
#' @param random Character or vector of characters that indicates the random effect factors.
#' @param modelname The name of the model file; optional.
#' @param wd The working directory in which OpenBUGS will work (i.e., generate the model files and chain information).
#' @return A model \sQuote{.txt} file is generated in the specified working directory. The function also returns a list of values:
#'  \describe{
#'   \item{init1,init2}{Default initial values for MCMC two chain procedure.}
#'   \item{vars}{A list of variables that are included in the estimation.}
#'   \item{nodes_sample}{a list of characters that specify the nodes to be monitored.}
#' }
#' 
#' @export
#' 
#' @seealso \code{\link{qrBugs}}
#' @examples
#'\dontrun{
#' # Need write access in the working directory before executing the code.
#' # No random component
#' bugsModel(y ~ x1 | x2, 't2','t2', random = NULL)
#' # Random component as subject ID
#' bugsModel(y ~ x1 | x2, 't2','t2', random = 'ID')
#'}
#'
bugsModel <- function(formula, fd, sd, random = NULL, modelname = "bugmodel", 
                      wd = getwd()) {
  
  fd <- tolower(fd)
  sd <- tolower(sd)
  
  # ********************** Prepare some default model parts
  
  modelpart01 <- "\n  model {\n  pi <- 3.14159265358979\n  for (i in 1:N) {  \n  dummy[i] <- 0\n  dummy[i] ~ dloglik(logLike[i])  \n  "
  
  # Likelihood part
  likelihood <- bugsLikelihood(fd, sd)
  
  # Initialize default sampling nodes - intercept only
  nodes_sample <- c("b_0", "d_0")
  
  # Initialize default initial values
  init1 <- list(0.1, 0.1)  # prepare initial values
  
  
  # ********************** Prepare model terms for composing the model main part
  formula <- as.Formula(formula)
  lxterm <- labels(terms(formula, rhs = 1L))  #extract mean predictors
  
  if (length(attr(formula, "rhs")) > 1) {
    # If dispersion model is specified and the supplied formula has the form y ~
    # X1+...+Xn | X1+...+Xn
    dxterm <- labels(terms(formula, rhs = 2L))  #extract dispersion predictors
  } else {
    # else dispersion model is not specified and the supplied formula has the form y ~
    # X1+...+Xn
    dxterm <- list()
  }
  
  
  # ********************** Start to compose model parts Mean/location model
  lm <- "mu[i] <- b_0"  # location model part -intercept
  lp_prior <- "b_0 ~dnorm(0.0, 1.0E-6)"  # location model part -intercept prior
  
  if (length(lxterm) != 0) {
    for (i in 1:length(lxterm)) {
      # add the location model terms to model part, priors, initial values and nodes
      # for sampling
      lm <- paste(lm, " + b_", gsub(":", "_", x = lxterm[i]), "*", gsub(":", 
        "[i]*", x = lxterm[i]), "[i]", sep = "")
      lp_prior <- paste(lp_prior, "\nb_", gsub(":", "_", x = lxterm[i]), "~dnorm(0.0, 1.0E-6)", 
        sep = "")
      init1 <- c(init1, 0.1)
      nodes_sample <- c(nodes_sample, paste("b_", gsub(":", "_", x = lxterm[i]), 
        sep = ""))
    }
  }
  
  # Dispersion model
  dm <- "log(sigma[i]) <- d_0"
  dp_prior <- "d_0 ~ dnorm(0.0, 1.0E-6)"
  
  if (length(dxterm) != 0) {
    for (i in 1:length(dxterm)) {
      # add the dispersion model terms to model part, priors, initial values and nodes
      # for sampling
      dm <- paste(dm, " + d_", gsub(":", "_", x = dxterm[i]), "*", gsub(":", 
        "[i]*", x = dxterm[i]), "[i]", sep = "")
      dp_prior <- paste(dp_prior, "\nd_", gsub(":", "_", x = dxterm[i]), "~dnorm(0.0, 1.0E-6)", 
        sep = "")
      init1 <- c(init1, 0.1)
      nodes_sample <- c(nodes_sample, paste("d_", gsub(":", "_", x = dxterm[i]), 
        sep = ""))
    }
  }
  
  
  # For random intercept part (currently only take the random effect for the mean
  # submodel)
  if (!is.null(random)) {
    
    formulas_random <- NULL
    
    for (i in 1:length(random)) {
      # add the random term to location submodel lm ~ .+u_XX[XX[i]]; XX can be a random
      # factor such as subject ID
      
      lm <- paste(lm, " + u_", random[i], "[", random[i], "[i]]", sep = "")
      
      
      # compose the random intercept estimation part for (j in 1:J) u_XX ~ dnorm(0,
      # tau_XX)
      formulas_random <- paste(formulas_random, "\n for (j in 1:J_", random[i], 
        ")\n {\n u_", random[i], "[j] ~ dnorm(0,tau_", random[i], ") \n", 
        sep = "")
      
      # The prior for tau tau_XX <- exp(2*lsu_XX); lsu_XX ~ dnorm(....)
      lp_prior <- paste(lp_prior, "\ntau_", random[i], " <- exp(2*lsu_", 
        random[i], ")\n lsu_", random[i], "~ dnorm(0.5, 1.0E-6)\n", sep = "")
      
      # Assume the initial value for lsu for each random term is 0.5
      init1 <- c(init1, 0.5)
      
      nodes_sample <- c(nodes_sample, paste("lsu_", random[i], sep = ""))
    }
    
    # Put mean and dispersion models/priors together
    formulas <- paste(lm, "\n", dm, "\n", "}\n", formulas_random, sep = "")
  }
  
  # Put mean and dispersion models/priors together if there is no random component
  if (is.null(random)) 
    formulas <- paste(lm, "\n", dm, "\n", sep = "")
  
  
  # Put priors together
  Priors <- paste("}\n", lp_prior, "\n", dp_prior, "\n", sep = "")
  
  # Put model front matter, likelihood, formulas, and priors together
  model <- paste(modelpart01, likelihood, "\n", formulas, Priors, "\n}", sep = "")
  
  # Generate the model file
  cat(model, file = paste(wd, "/", modelname, ".txt", sep = ""))
  
  # Produce default initial value list
  names(init1) <- nodes_sample
  init2 <- init1
  
  
  list(init1 = init1, init2 = init2, nodes_sample = nodes_sample, 
       vars = attr(terms(formula), 
    "variables"))
  
} 
