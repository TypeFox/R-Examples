#-----------------------------------------------------------------------------#
#' Estimate Causal Effects in presence of interference  
#'
#' @param formula The formula used to define the causal model.  Has a minimum 
#' of 4 parts, separated by \code{|} and \code{~} in a specific structure: 
#' \code{outcome | exposure ~ propensity covariates | group}. The order matters, 
#' and the pipes split the data frame into corresponding pieces. The part 
#' separated by \code{~} is passed to the chosen \code{model_method} used to 
#' estimate or fix propensity parameters. 
#' @param propensity_integrand A function, which may be created by the user, 
#' used to compute the IP weights. This defaults to \code{logit_integrand}, 
#' which calculates the product of inverse logits for individuals in a group: 
#' \eqn{\prod_{j = 1}^{n_i} \{r \times h_{ij}(b_i)^{A_{ij}}\}\{1 - r \times 
#' h_{ij}(b_i)^{1 - A_{ij}} \} f_b(b_i; \theta_s)}{prod(r * plogis(X * fixef + b)^A * 
#' (1 - r * plogis(X * fixef+ b))^(1 - A)) * 
#' dnorm(sd = sqrt(ranef))} where \deqn{h_{ij}(b_i) = logit^{-1}
#' (\mathbf{X}_{ij}\theta_a + b_i)} and \eqn{b_i} is a group-level random effect, 
#' \eqn{f_b} is a \eqn{N(0, \theta_s)} density, and \eqn{r} is a known 
#' randomization probability which may be useful if a participation vector is 
#' included in the \code{formula}. If no random effect was included in the 
#' \code{formula}, \code{logit_integrand} essentially ignores the random effect 
#' and \eqn{f_b(b_i, \theta_s)} integrates to 1. See details for arguments that 
#' can be passed to \code{logit_integrand}
#' @param loglihood_integrand A function, which may be created by the user, that
#' defines the log likelihood of the logit model used for \code{robust} variance
#' estimation. Generally, this will be the same function as 
#' \code{propensity_integrand}. Indeed, this is the default.
#' @param allocations a vector of values in [0, 1]. Increasing the number of 
#' elements of the allocation vector greatly increases computation time; however, 
#' a larger number of allocations will make plots look nicer. A minimum of two 
#' allocations is required.
#' @param data the analysis data frame. This must include all the variables 
#' defined in the \code{formula}.
#' @param model_method the method used to estimate or set the propensity model 
#' parameters. Must be one of \code{'glm'}, \code{'glmer'}, or \code{'oracle'}. 
#' Defaults to \code{'glmer'}. For a fixed effects only model use \code{'glm'}, 
#' and to include random effects use\code{'glmer'}. \code{logit_integrand} 
#' only supports a single random effect for the grouping variable, so if more 
#' random effects are included in the model, different \code{propensity_integrand}
#' and \code{loglihood_integrand} functions should be defined. When the propensity 
#' parameters are known (as in simulations) or if estimating parameters by 
#' other methods, use the \code{'oracle'} option. See \code{model_options} for
#' details on how to pass the oracle parameters. 
#' @param model_options a list of options passed to the function in 
#' \code{model_method}. Defaults to \code{list(family = binomial(link = 'logit'))}. 
#' When \code{model_method = 'oracle'}, the list must have two elements 
#' \code{fixed.effects} and \code{random.effects}. If the model did not include 
#' random effects, set \code{random.effects = NULL}.
#' @param causal_estimation_method currently only supports \code{'ipw'}.
#' @param causal_estimation_options A list with two slots. (1) \code{variance_estimation} is 
#' either \code{'naive'} or \code{'robust'}. See details. Defaults to \code{'robust'}. (2) 
#' is \code{set_NA_to_0}. Defaults to \code{TRUE}. When, for example, group sizes 
#' reach over 1000, the product terms of the propensity diminish to zero. 
#' This may result in \code{NaN} values for the weights or loglihood. This option
#' sets such cases to zero.
#' @param conf.level level for confidence intervals. Defaults to \code{0.95}.
#' @param rescale.factor a scalar multiplication factor by which to rescale outcomes
#' and effects. Defaults to \code{1}.
#' @param ... Used to pass additional arguments to internal functions such as 
#' \code{numDeriv::grad()} or \code{integrate()}. Additionally, arguments can be 
#' passed to the \code{propensity_integrand} and \code{loglihood_integrand} functions.
#'
#' @details The following formula includes a random effect for the group: \code{outcome | 
#' exposure ~ propensity covariates + (1|group) | group}. In this instance, the 
#' group variable appears twice. If the study design includes a "participation" 
#' variable, this is easily added to the formula: \code{outcome | exposure | 
#' participation ~ propensity covariates | group}. 
#' 
#' \code{logit_integrand} has two options that can be passed via the \code{...} 
#' argument: 
#' \itemize{
#'  \item \code{randomization}: a scalar. This is the \eqn{r} in the formula just 
#'  above. It defaults to 1 in the case that a \code{participation} vector is not 
#'  included. The vaccine study example demonstrates use of this argument.
#'  \item \code{integrate_allocation}: \code{TRUE/FALSE}. When group sizes grow 
#'  large (over 1000), the product term of \code{logit_integrand} tends quickly to 0. 
#'  When set to \code{TRUE}, the IP weights tend less quickly to 0. 
#'  Defaults to \code{FALSE}.
#' }
#' 
#' If the true propensity model is known (e.g. in simulations) use 
#' \code{variance_estimatation = 'naive'}; otherwise, use the default 
#' \code{variance_estimatation = 'robust'}. Refer to the web appendix of
#' \href{http://dx.doi.org/10.1111/biom.12184}{Heydrich-Perez et al. (2014)} 
#' for complete details.
#' 
#' 
#' @return Returns a list of overall and group-level IPW point estimates, overall 
#' and group-level IPW point estimates (using the weight derivatives), derivatives
#' of the loglihood, the computed weight matrix, the computed 
#' weight derivative array, and a summary.
#' 
#' @export
#' 
#-----------------------------------------------------------------------------#

interference <- function(formula,
                         propensity_integrand = 'logit_integrand',
                         loglihood_integrand = propensity_integrand,
                         allocations,
                         data,
                         model_method = "glmer",
                         model_options = list(family = binomial(link = 'logit')),
                         causal_estimation_method = 'ipw',
                         causal_estimation_options = list(set_NA_to_0 = TRUE, 
                                                          variance_estimation = 'robust'),
                         conf.level = 0.95,
                         rescale.factor = 1,   
                         ...)
{
  ## Necessary bits ##
  dots <- list(...)
  integrandFUN    <- match.fun(propensity_integrand)
  likelihoodFUN   <- match.fun(loglihood_integrand)
  oracle          <- model_method == 'oracle'
  cformula        <- Formula::Formula(formula)
  len_lhs         <- length(cformula)[1]
  len_rhs         <- length(cformula)[2]

  ## For the sake of consistency downstream, reorder data frame by groups ##
  group_var <- attr(terms(cformula, lhs = 0, rhs = len_rhs), 'term.labels')
  data <- data[order(data[ , group_var]), ]
  
  ## Parse out the formula into necessary pieces ##
  Y <- Formula::model.part(cformula, data = data, lhs = 1, drop = TRUE)
  A <- Formula::model.part(cformula, data = data, lhs = 2, drop = TRUE)
  G <- Formula::model.part(cformula, data = data, rhs = len_rhs, drop = TRUE)
  
  # Used when there is 'participation' variable
  if(len_lhs > 2){
    B <- Formula::model.part(cformula, data = data, lhs = len_lhs, drop = TRUE)
  } else {
    B <- A
  }
  
  propensity_formula <- formula(terms(cformula, lhs = len_lhs, rhs = -2))
  random.count <- length(lme4::findbars(propensity_formula))

  trt_lvls     <- sort(unique(A))
  N            <- length(unique(G))
  k            <- length(allocations)
  l            <- length(trt_lvls)

  ## Warnings ##
  if(model_method == 'glm' & random.count > 0 ){
    stop('propensity_formula appears to include a random effect when "glm" was chosen \n 
         for parameter estimation. Set model_method to "glmer" to include a random effect')
  }
  
  if(propensity_integrand == "logit_integrand" & random.count > 1){
    stop('Logit integrand is designed to handle only 1 random effect.')
  }
  
  if(min(allocations) < 0 | max(allocations) > 1){
    stop('Allocations must be between 0 and 1 (inclusive)')
  }
  
  if(length(allocations) < 2){
    warning('At least 2 allocations must be specified in order to estimate indirect effects')
  }
  
  if(N == 1){
    stop('The group variable must have at least 2 groups (more is better).')
  }
  
  #### Compute Parameter Estimates ####

  estimation_args <- append(list(formula = propensity_formula, data = data), 
                            model_options)
  
  if(model_method == "glmer"){
    propensity_model <- do.call(lme4::glmer, args = estimation_args)
    fixed.effects  <- lme4::getME(propensity_model, 'fixef')
    random.effects <- lme4::getME(propensity_model, 'theta')
    X <- lme4::getME(propensity_model, "X")
  } else if(model_method == "glm"){
    propensity_model <- do.call("glm", args = estimation_args)
    fixed.effects  <- coef(propensity_model)
    random.effects <- NULL
    X <- model.matrix(propensity_model)
  } else if(model_method == "oracle"){
    fixed.effects  <- model_options[[1]]
    random.effects <- model_options[[2]]
    X <- model.matrix(propensity_formula, data)
    
    if(length(fixed.effects) != ncol(X)){
      stop('oracle fixed effects vector must have length of # of columns of X')
    }
  }
  
  #### Compute Effect Estimates ####
  out <- list()
  grid <- effect_grid(allocations = allocations, treatments  = trt_lvls)
  
  # Will have other causal estimation types in the future
  if('ipw' %in% causal_estimation_method)
  {
    ipw_args <- append(append(dots, causal_estimation_options),
                       list(propensity_integrand = integrandFUN, 
                            loglihood_integrand  = likelihoodFUN,
                            allocations          = allocations,
                            fixed.effects        = fixed.effects, 
                            random.effects       = random.effects,
                            Y = Y, X = X, A = A, B = B, G = G))
  
    ipw <- do.call(ipw_interference, args = ipw_args)
    out <- append(out, ipw)
    
    print('Computing effect estimates...')
    
    estimate_args <- list(obj = ipw,
                          variance_estimation = causal_estimation_options$variance_estimation,
                          alpha1      = grid$alpha1,
                          trt.lvl1    = grid$trt1,
                          alpha2      = grid$alpha2,
                          trt.lvl2    = grid$trt2,
                          marginal    = grid$marginal,
                          effect_type = grid$effect_type,
                          rescale.factor = rescale.factor,
                          conf.level = conf.level,
                          print = FALSE)
    
    est <- do.call(ipw_effect_calc, args = estimate_args)
    
    # 2/21/15: ipw_effect_calc returns a data.frame of lists, but these should 
    # be numeric columns. Here's a quick fix.
    est <- apply(t(est), 2, as.numeric)
    
    out$estimates <- cbind(grid, est)
  }

  #### Summary ####
  
  # count missing weights by allocation
  weights_na <- apply(out$weights, 2, function(x) sum(is.na(x)))
  
  if(!oracle){
    out$models <- list(propensity_model = propensity_model)
  } else {
    out$models <- list(propensity_model = model_options) 
  }
  
  out$summary <- list(formula      = deparse(formula),
                      oracle       = oracle,
                      conf.level   = conf.level,
                      ngroups      = N, 
                      nallocations = k,
                      npredictors  = length(fixed.effects),
                      ntreatments  = l,
                      allocations  = allocations,
                      treatments   = trt_lvls,
                      weights_na_count  = weights_na,
                      variance_estimation = causal_estimation_options$variance_estimation)
  
  class(out) <- "interference"
  
  print('Interference complete')
  return(out)
}