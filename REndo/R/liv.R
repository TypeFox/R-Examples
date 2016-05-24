#'@title  Fitting Linear Models with one Endogenous Regressor using Latent Instrumental Variables
#'@aliases liv
# Description
#'@description  Fits linear models with one endogenous regressor and no additional explanatory variables using the latent instrumental variable approach
#'presented in Ebbes,P., Wedel,M.,  B\"{o}ckenholt, U., and Steerneman, A. G. M. (2005). This is a statistical technique to address the endogeneity problem where no external instrumental
#'variables are needed. The important assumption of the model is that the latent variables are discrete with at least two groups with different means and
#'the structural error is normally distributed.
#
# Arguments
#'@param    formula  an object of type 'formula': a symbolic description of the model to be fitted. Example \code{var1 ~ var2}, where \code{var1} is a vector
#' containing the dependent variable, while \code{var2} is a vector containing the endogenous variable.
#'@param    data  optional data frame or list containing the variables of the model.
#'@param    param  a vector of initial values for the parameters of the model to be supplied to the optimization algorithm. In any model there are eight parameters.
#'The first parameter is the intercept, then the coefficient of the endogenous variable followed by the means of the two groups of the latent IV (they need to be different, otherwise model is not identified),
#'then the next three parameters are for the variance-covariance matrix. The last parameter is the probability of being in group 1. When not provided, 
#'initial paramameters values are set equal to the OLS coefficients, the two group means are set to be equal to \code{mean(P)} and \code{ mean(P) + sd(P)}, the
#'variance-covariance matrix has all elements equal to 1 while \code{probG1} is set to equal 0.5.
#'
#'@details     Let's consider the model:
#' \deqn{Y_{t} = \beta_{0} + \alpha P_{t} + \epsilon_{t}}{Y_t = b0 + a * P_t + eps_t}
#' \deqn{P_{t}=\pi^{'}Z_{t} + \nu_{t}}{P_t = pi * Z_t + nu_t}
#' where \code{t = 1,..,T} indexes either time or cross-sectional units, \eqn{Y_{t}}{Y_t} is the dependent variable, \eqn{P_{t}}{P_t} is a \code{k x 1} continuous, endogenous regressor, 
#' \eqn{\epsilon_{t}}{epsilon_t} is a structural error term with mean zero and \eqn{E(\epsilon^{2})=\sigma^{2}_{\epsilon}}{E(eps^2) = sigma_eps^2}, \eqn{\alpha}{a} and \eqn{\beta}{b0}
#' are model parameters. \eqn{Z_{t}}{Z_t} is a \code{l x 1} vector of instruments, and \eqn{\nu_{t}}{nu} is a random error with mean zero and \eqn{E(\nu^{2}) = \sigma^{2}_{\nu}}{E(nu^2) = sigma_nu^2}. 
#' The endogeneity problem arises from the correlation of \code{P} and \eqn{\epsilon_{t}}{eps} through \eqn{E(\epsilon\nu) = \sigma_{\epsilon\nu}}{E(eps * nu) = sigma_0^2}.
#' 
#' LIV  considers \eqn{Z_{t}^{'}}{Z_t} to be a latent, discrete, exogenous variable with an unknown number of groups \code{m} and \eqn{\pi}{pi} is a vector of group means. 
#' It is assumed that \code{Z} is independent of the error terms \eqn{\epsilon}{eps} and \eqn{\nu}{nu} and that it has at least two groups with different means. 
#' The structural and random errors are considered normally distributed with mean zero and variance-covariance matrix \eqn{\Sigma}{Sigma}: 
#' \deqn{\Sigma = \left(
#' \begin{array}{ccc}
#' \sigma_{\epsilon}^{2} & \sigma_{\epsilon\nu}\\
#' \sigma_{\epsilon\nu} & \sigma_{\nu}^{2}
#' \end{array}\right)}{Sigma = (sigma_eps^2,       sigma_0^2 ; sigma_0^2,       sigma_nu^2)}
#' The identification of the model lies in the assumption of the non-normality of \eqn{P_{t}}{P}, the discreteness of the unobserved instruments and the existence of
#' at least two groups with different means. 
#' 
#' The method has been programmed such that the latent variable has two groups. Ebbes et al.(2005) show in a Monte Carlo experiement that
#' even if the true number of the categories of the instrument is larger than two, LIV estimates are approximately consistent. Besides, overfitting in terms
#' of the number of groups/categories reduces the degrees of freedom and leads to efficiency loss. When provided by the user, the initial parameter values
#' for the two group means have to be different, otherwise the model is not identified. For a model with additonal explanatory variables a Bayesian approach is needed, since
#' in a frequentist approach identification issues appear. The optimization algorithm used is BFGS.
#'
#Return Value
#'@return    Returns the optimal values of the parameters as computed by maximum likelihood using BFGS algorithm. 
#'\item{coefficients}{returns the value of the parameters for the intercept and the endogenous regressor as computed with maximum likelihood.}
#'\item{means}{returns the value of the parameters for the means of the two categories/groups of the latent instrumental variable.}
#'\item{sigma}{returns the variance-covariance matrix sigma, where on the main diagonal are the variances of the structural error and that of
#'the endogenous regressor and the off-diagonal terms are equal to the covariance between the errors.}
#'\item{probG1}{returns the probability of being in group one. Since the model assumes that the latent instrumental variable has two groups,
#'\code{1-probG1} gives the probability of group 2.}
#'\item{value}{the value of the log-likelihood function corresponding to the optimal parameters.}
#'\item{convcode}{an integer code, the same as the output returned by \code{optimx}. 0 indicates successful completion. A possible error code is 1 which indicates that the iteration
#'limit maxit had been reached.}
#'\item{hessian}{a symmetric matrix giving an estimate of the Hessian at the solution found.}
#'@keywords endogenousdata
#'@keywords latent
#'@keywords instruments
#'@author The implementation of the model formula by Raluca Gui based on the paper of Ebbes et al. (2005).
#'@references   Ebbes, P., Wedel,M., B\"{o}ckenholt, U., and Steerneman, A. G. M. (2005). 'Solving and Testing for Regressor-Error
#'(in)Dependence When no Instrumental Variables are Available: With New Evidence for the Effect of Education on Income'. 
#'\emph{Quantitative Marketing and Economics},
#' \bold{3}:365--392.
#' @examples
#' # load data
#' data(dataLIV)
#' y <- dataLIV$y
#' P <- dataLIV$P
#' # function call without any initial parameter values 
#' l  <- liv(y ~ P)
#' summary(l)
#' # function call with initial parameter values given by the user
#' l1 <- liv(y ~ P, c(2.9,-0.85,0,0.1,1,1,1,0.5))
#' summary(l1)
# make availble to the package users
#'@export
liv <- function(formula, param=NULL, data=NULL){

 if( ncol(stats::get_all_vars(formula)) != 2 )
    stop("A wrong number of parameters were passed in the formula. No exogenous variables are admitted.")

  mf <- stats::model.frame(formula = formula, data = data)

  
  # if user parameters are not defined, provide initial param. values
  # coefficients are the OLS coefficients
  # the two group means = mean(P), and mean(P) + sd(P)
  # next three parameters = 1
  # probG1 = 0.5
  
  y <- mf[,1]
  P <- mf[,ncol(mf)]
  
  if (is.null(param)) {

    param1 <- stats::coefficients(stats::lm(mf[,1]~mf[,2]))[1]
    param2 <- stats::coefficients(stats::lm(mf[,1]~mf[,2]))[2]
    param3 <- mean(mf[,2])
    param4 <- mean(mf[,2]) + stats::sd(mf[,2])
    param5 <- param6 <- param7 <- 1
    param8 <- 0.5
    param <- as.double(c(param1,param2,param3,param4,param5,param6,param7,param8))
  }
  
  b <- optimx::optimx( par=param,fn=logL, y=mf[,1],P=mf[,2],
                       method="BFGS",hessian=T,
                       control=list(trace=0))
  
  obj <- methods::new("liv")
  obj@formula <- formula
  
  obj@coefficients <- c(b$p1, b$p2)      # coefficients
  
  obj@groupMeans <-  c(b$p3, b$p4)     # means of the 2 groups of the latent IV
  
  # variance-covariance matrix of errors
  obj@sigma <- matrix(c(b$p5^2,b$p5*b$p6,b$p5*b$p6,b$p6^2+b$p7^2),2,2)
  
  obj@probG1 <- b$p8     # probability of group 1
  # check for the probability to be less than 1 - if not, give warning, change initial parameter values
  if (obj@probG1 > 1) warning("Probability of Group 1 greater than 0. Check initial parameter values")
  
  obj@initValues <- param 
  obj@value <- b$value          # the value of the likelihood function corresponding to param
  obj@convCode <- as.integer(b$convcode)    # message whether if converged
  
  hess <- attr(b,"details")[,"nhatend"]         # hessian matrix
  std_par <- suppressWarnings(sqrt(diag(solve(do.call(rbind,hess)))))
  obj@hessian <- hess[[1]]
  
  
  obj@seCoefficients <- std_par[1:2]
  if (obj@seCoefficients[1] =="NaN") warning("Coefficient's standard errors unable to be computed. Check initial parameter values")
  
  obj@seMeans <- std_par[3:4]
  if (obj@seMeans[1] =="NaN") warning("Group means' standard errors unable to be computed. Check initial parameter values")
  
  obj@seProbG1 <- std_par[8]
  
  return(obj)
}
