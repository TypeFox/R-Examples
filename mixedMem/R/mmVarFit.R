#' Fit  Mixed Membership models using variational EM
#' 
#' \code{mmVarFit} is the primary function of the \code{mixedMem} package. The function fits parameters \eqn{\phi} and \eqn{\delta} for the variational
#'  distribution of latent variables as well as psuedo-MLE estimates for 
#'  the population parameters \eqn{\alpha} and \eqn{\theta}. See documentation for
#'  \code{mixedMemModel} or the package vignette, "Fitting Mixed Membership Models Using \code{mixedMem}" for a more detailed description of 
#'  the variables and notation in a mixed membership model. 
#' 
#'  
#' \code{mmVarFit} selects psuedo-MLE estimates for \eqn{\alpha} and \eqn{\theta} and approximates
#' the posterior distribution for the latent variables through a mean field variational approach.
#' The variational lower bound on the log-likelihood at the data given model parameters is given by Jensen's inequality:
#' 
#' \eqn{E_Q{\log[p(X,Z, \Lambda)]} - E_Q{\log[Q(Z, \Lambda|\phi, \delta)]} \le \log P(obs |\alpha, \theta)} where
#' 
#' \eqn{Q = \prod_i [Dirichlet(\lambda_i|\phi_i) \prod_j^J \prod_r^{R_j} \prod_n^{N_{i,j,r}}Multinomial(Z_{i,j,r,n}|\delta_{i,j,r,n})]}
#' 
#' The left hand side of the above equation is often referred to as the Evidence Lower Bound (ELBO). The global parameters
#' \eqn{\alpha} and \eqn{\theta} are selected to maximize this bound as a psuedo-MLE procedure. In addition, it can be shown that
#' maximizing the ELBO with respect to the variational parameters \eqn{\phi} and \eqn{\delta} is equivalent to   
#' minimizing the KL divergence between the tractable variational distribution and the true posterior.
#'
#' The method uses a variational EM approach. The E-step considers the global parameters \eqn{\alpha} and \eqn{\theta} fixed
#' and picks appropriate variational parameters \eqn{\phi} and \eqn{\delta} to minimize the KL divergence between the 
#' true posterior and the variational distribution (or maximize the ELBO). The M-step considers
#' the variational parameters  fixed and global parameters \eqn{\theta} and \eqn{\alpha} are selected which maximize the lower bound 
#' on the log-likelihood (ELBO).
#'   
#' @references
#' Beal, Matthew James. Variational algorithms for approximate Bayesian inference. Diss. University of London, 2003.
#'
#' @param model a \code{mixedMemModel} object created by the \code{mixedMemModel} constructor.
#' @param printStatus an integer 0 or 1. When \code{printStatus} is 1 \code{mmVarFit} will print status updates, when
#' \code{printStatus} is 0 \code{mmVarFit} will not print status updates.
#' @param printMod a positive integer which specifies how often to print status updates. The 
#' status will be printed at each step which is a multiple of \code{printMod}. The default value is 1.
#' @param stepType an integer from 0-3 which specifies what steps to carry out in estimation. 0 performs a single E-Step; this can be used
#' to find the held out ELBO. 1 performs iterates between E-steps and M-steps for estimating theta but keeps alpha constant.
#'  2 iterates between  E-steps and M-steps for estimating alpha, but keeps theta constant. 3 completes full
#' E and M steps. The default is 3.
#' @param maxTotalIter the maximum total steps before termination. A full E and M step together count as 1 step.
#'  If this maximum is ever achieved, a warning message will be printed at convergence. The default is 500.
#' @param maxEIter the maximum iterations for each E-Step. If this maximum is ever achieved, a warning 
#' message will be printed at convergence. The default is 1000.
#' @param maxAlphaIter the maximum iterations when fitting alpha. If this maximum is ever achieved, 
#' a warning message will be printed at convergence. The default is 200.
#' @param maxThetaIter the maximum iterations when fitting theta. If this maximum is ever achieved, 
#' a warning message will be printed at convergence. The default is 1000.
#' @param maxLSIter the maximum backtracking iterations in the line search for updating alpha and theta for rank data. 
#' @param elboTol the convergence criteria for the EM Algorithim. When the relative increase in the ELBO is less than the convergence critiera,
#' the algorithm converges. The default is 400.
#' @param alphaTol the convergence criteria for updates to alpha. When the relative increase in the ELBO is less than the convergence critiera,
#' the update to alpha converges. The defualt is 1e-6.
#' @param thetaTol the convergence criteria for updates to theta. When the relative increase in the ELBO is less than the convergence critiera,
#' the update to theta converges. The default is 1e-10.
#' @param aNaught the first step size in the backtracking line search used to update alpha and theta for rank data. The default is 1.
#' @param tau the backtracking factor in the backtracking line search used to update alpha and theta for rank data. The default is .899.
#' @param bMax the number of times the penalization factor is increased in the interior point method 
#' for fitting theta for rank data. The default is 3.
#' @param bNaught the initial penalization factor in the interior point method for fitting theta for rank data. 
#' \code{bNaught} must be positive and the default value is 1000.
#' @param bMult a positive number by which bNaught is multiplied by in each iteration of the interior point method 
#' for fitting theta for rank data. The default is 1000.
#' @param vCutoff a positive integer cutoff for Vj at which a gradient ascent method is used instead of the Newton Raphson interior point method. This is used to avoid inverting
#' a large matrix. The default is 13.
#' @param holdConst a vector of integers specifying groups to be held fixed during the estimation procedure. The estimation algorithim will hold the
#'  theta parameters of these specific groups constant, but update all other parameters. The group numbers range from 0 to K - 1. To estimate all groups, use the default value of c(-1). 
#' @return \code{mmVarFit} returns a \code{mixedMemModel} object containing updated variational parameters and global parameters.
#' @seealso mixedMemModel
#' @examples
#' ## The following example generates multivariate observations (both variables are multinomial)
#' ## from a mixed membership model. The observed data is then used to fit a 
#' ## mixed membership model using mmVarFit
#' 
#' ## Generate Data 
#' Total <- 30 #30 Individuals
#' J <- 2 # 2 variables
#' dist <- rep("multinomial",J) # both variables are multinomial
#' Rj <- rep(100,J) # 100 repititions for each variable
#' # Nijr will always be 1 for multinomials and bernoulli's
#' Nijr <- array(1, dim = c(Total, J, max(Rj))) 
#' K <- 4 # 4 sub-populations
#' alpha <- rep(.5, K) # hyperparameter for dirichlet distribution
#' Vj <- rep(5, J) # each multinomial has 5 categories
#' theta <- array(0, dim = c(J, K, max(Vj)))
#' theta[1,,] <- gtools::rdirichlet(K, rep(.3, 5))
#' theta[2,,] <- gtools::rdirichlet(K, rep(.3, 5))
#' obs <- rmixedMem(Total, J, Rj, Nijr, K, Vj, dist, theta, alpha)$obs
#' 
#' ## Initialize a mixedMemModel object
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj, Nijr= Nijr,
#'  K = K, Vj = Vj,dist = dist, obs = obs,
#'   alpha = alpha, theta = theta+0)
#' 
#' ## Fit the mixed membership model
#' out <-mmVarFit(test_model)
#' @export
mmVarFit = function(model, printStatus = 1,
                    printMod = 1, stepType = 3,
                    maxTotalIter = 500, maxEIter = 1000,
                    maxAlphaIter = 200, maxThetaIter = 1000,
                    maxLSIter = 400, elboTol = 1e-6, alphaTol = 1e-6,
                    thetaTol = 1e-10, aNaught = 1.0, tau = .899,
                    bMax = 3, bNaught = 1000.0, bMult = 1000.0, vCutoff = 13, holdConst = c(-1)) {
  
  checkModel(model) # R function which checks inputs
  print("Model Check: Ok!")
  output <- model
  names(output) = c("Total", "J", "Rj", "Nijr", "K", "Vj", "alpha",
                    "theta", "phi", "delta", "dist" ,"obs",
                    "fixedObs", "P", "beta")[1:length(model)]
  if(is.null(model$fixedObs)){
    print("<== Beginning Model Fit! ==>")
    ret <- varInfInputC(model, printStatus, printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter,
                 maxThetaIter, maxLSIter, elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, 
                 bMult, vCutoff, holdConst) # R wrapper function
    output$alpha <- ret$alpha
    output$theta <- ret$theta
    output$phi <- ret$phi
    output$delta <- ret$delta
    
  } else {
    print("<== Beginning Extended Model Fit! ==>")
    ret <- varInfInputExtC(model, printStatus, printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter,
                 maxThetaIter, maxLSIter, elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, 
                 bMult, vCutoff, holdConst) # R wrapper function
    output$alpha <- ret$alpha
    output$theta <- ret$theta
    output$phi <- ret$phi
    output$delta <- ret$delta
    output$P <- ret$P
    output$beta <- ret$beta
  }

  return(output)
}



#' Compute a lower bound on the log-likelihood (ELBO)
#' 
#' \code{computeELBO} computes the variational lower bound on the log-likelihood, also called the ELBO, for a mixed membership model.
#'
#' The lower bound (ELBO) is the objective function in the variational EM algorithim. 
#' It is a function of the latent variables (\eqn{\phi} and \eqn{\delta}) and the parameters (\eqn{\alpha} and \eqn{\theta})
#' that be derived from Jensen's inequality:
#'      
#' \eqn{E_Q{\log[p(X,Z, \Lambda)]} - E_Q{\log[Q(Z, \Lambda|\phi, \delta)]} \le\log P(obs |\alpha, \theta) }
#'    
#' @param model a \code{mixedMemModel} object created by the \code{mixedMemModel} constructor.
#' @return \code{computeELBO} returns the lower bound on the log-likelihood, a real number.
#' @export
computeELBO = function(model)
{
  checkModel(model)
  if(is.null(model$P)){
    return(computeElboC(model))
  } else {
    return(computeElboExtC(model))
  }
}