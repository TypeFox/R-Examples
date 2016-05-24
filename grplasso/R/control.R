setClass("grpl.control",
         representation = representation(
           save.x       = "logical",
           save.y       = "logical",
           update.hess  = "character",
           update.every = "numeric",
           inner.loops  = "numeric",
           line.search  = "logical",
           max.iter     = "numeric",
           tol          = "numeric",
           lower        = "numeric",
           upper        = "numeric",
           beta         = "numeric",
           sigma        = "numeric",
           trace        = "numeric"),

         prototype = list(
           save.x       = FALSE,
           save.y       = TRUE,
           update.hess  = "lambda",
           update.every = 3,
           inner.loops  = 10,
           line.search  = TRUE,
           max.iter     = 500,
           tol          = 5 * 10^-8,
           lower        = 10^-2,
           upper        = 10^9,
           beta         = 0.5,
           sigma        = 0.1,
           trace        = 1),
           
         validity = function(object){
           if(ceiling(object@update.every) != floor(object@update.every) |
              object@update.every <= 0)
             return("update.every has to be a natural number")

           if(ceiling(object@inner.loops) != floor(object@inner.loops) |
              object@inner.loops < 0)
             return("inner.loops has to be a natural number or 0")

           if(ceiling(object@max.iter) != floor(object@max.iter) |
              object@max.iter <= 0)
             return("inner.loops has to be a natural number or greater than 0")

           if(object@beta <= 0 | object@beta >= 1)
             return("beta has to be in (0, 1)")
           
           if(object@sigma <= 0 | object@sigma >= 1)
             return("sigma has to be in (0, 1)")

           if(object@tol <= 0)
             return("tol has to be positive")

           if(object@lower > object@upper)
             return("lower <= upper has to hold")

           return(TRUE)
         }
)

grpl.control <- function(save.x = FALSE, save.y = TRUE,
                         update.hess = c("lambda", "always"),
                         update.every = 3, inner.loops = 10,
                         line.search = TRUE, max.iter = 500,
                         tol = 5 * 10^-8, lower = 10^-2, upper = Inf,
                         beta = 0.5, sigma = 0.1, trace = 1){
  
  ## Purpose: Options for the Group Lasso Algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## save.x: a logical indicating whether the design matrix should be saved.
  ## save.y: a logical indicating whether the response should be saved.
  ## update.hess: should the hessian be updated in each
  ##              iteration ("always")? update.hess = "lambda" will update
  ##              the Hessian once for each component of the penalty
  ##              parameter "lambda" based on the parameter estimates
  ##              corresponding to the previous value of the penalty
  ##              parameter.
  ## inner.loops: how many loops should be done (at maximum) when solving
  ##              only the active set (without considering the remaining
  ##              predictors)
  ## tol: convergence tolerance; the smaller the more precise, see
  ##      details below.
               
  ## lower: lower bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## upper: upper bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## beta: scaling factor beta < 1 of the Armijo line search.
  ## sigma: 0 < \sigma < 1 used in the Armijo line search.
  ## trace: integer. "0" omits any output,
  ##        "1" prints the current lambda value,
  ##        "2" prints the improvement in the objective function after each
  ##        sweep through all the parameter groups and additional
  ##        information.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:02

  
  update.hess <- match.arg(update.hess)

  RET <- new("grpl.control",
             save.x       = save.x,
             save.y       = save.y,
             update.hess  = update.hess,
             update.every = update.every,
             inner.loops  = inner.loops,
             line.search  = line.search,
             max.iter     = max.iter,
             tol          = tol,
             lower        = lower,
             upper        = upper,
             beta         = beta,
             sigma        = sigma,
             trace        = trace)
  RET
}
