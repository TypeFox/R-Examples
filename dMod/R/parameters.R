#' Parameter transformation
#' 
#' @param trafo Named character vector. Names correspond to the parameters being fed into
#' the model (the inner parameters). The elements of tafo are equations that express 
#' the inner parameters in terms of other parameters (the outer parameters)
#' @param parameters Character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(trafo)} the identity transformation is assumed.
#' @param compile Logical, compile the function (see \link{funC0})
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @seealso \link{Pi} for implicit parameter transformations and
#' \link{concatenation} for the concatenation of parameter transformations
#' @examples
#' \dontrun{
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' }
#' @export
P <- function(trafo, parameters=NULL, compile = FALSE) {
  
  # get outer parameters
  symbols <- getSymbols(trafo)
  
  if(is.null(parameters)) {
    parameters <- symbols 
  } else {
    identity <- parameters[which(!parameters%in%symbols)]
    names(identity) <- identity
    trafo <- c(trafo, identity)
  }
  
  
  # expression list for parameter and jacobian evaluation
  trafo.list <- lapply(trafo, function(myrel) parse(text=as.character(myrel)))
  jacobian <- unlist(lapply(parameters, function(var) {
    unlist(lapply(trafo.list, function(myexp) paste(deparse(D(myexp, as.character(var))), collapse="")))
  }))
  
  jacNames <- expand.grid.alt(names(trafo), parameters)
  jacNames <- paste(jacNames[,1], jacNames[,2], sep=".")
  
  dtrafo <- jacobian; names(dtrafo) <- jacNames
  
  PEval <- funC0(trafo, compile = compile)
  dPEval <- funC0(dtrafo, compile = compile)
  
  # the parameter transformation function to be returned
  p2p <- function(p, fixed=NULL, deriv = TRUE) {
    
    # Inherit from p
    dP <- attr(p, "deriv", exact = TRUE)
    
    # Evaluate transformation
    args <- c(as.list(p), as.list(fixed))
    pinner <- PEval(args)[1,]
    dpinner <- dPEval(args)[1,]
    
    # Construct output jacobian
    jac.vector <- rep(0, length(pinner)*length(p))
    names(jac.vector) <- outer(names(pinner), names(p), function(x, y) paste(x, y, sep = "."))
    
    names.intersect <- intersect(names(dpinner), names(jac.vector))
    jac.vector[names.intersect] <- as.numeric(dpinner[names.intersect])
    jac.matrix <- matrix(jac.vector, length(pinner), length(p), dimnames = list(names(pinner), names(p)))
    
    
    if(!is.null(dP)) jac.matrix <- jac.matrix%*%dP[colnames(jac.matrix),]
    if(deriv) attr(pinner, "deriv") <- jac.matrix
    
    return(pinner)
  }
  
  class(p2p) <- "par" 
  attr(p2p, "equations") <- trafo
  attr(p2p, "parameters") <- parameters
  
  
  return(p2p)
  
}



#' Parameter transformation (implicit)
#' 
#' @param trafo Named character vector defining the equations to be set to zero. 
#' Names correspond to dependent variables.
#' @param parameters Character vector, the independent variables.  
#' @param compile Logical, compile the function (see \link{funC0})
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @details Usually, the equations contain the dependent variables, the independent variables and 
#' other parameters. The argument \code{p} of \code{p2p} must provide values for the independent
#' variables and the parameters but ALSO FOR THE DEPENDENT VARIABLES. Those serve as initial guess
#' for the dependent variables. The dependent variables are then numerically computed by 
#' \link[rootSolve]{multiroot}. The Jacobian of the solution with respect to dependent variables
#' and parameters is computed by the implicit function theorem. The function \code{p2p} returns
#' all parameters as they are with corresponding 1-entries in the Jacobian.
#' #' @seealso \link{P} for explicit parameter transformations and
#' \link{concatenation} for the concatenation of parameter transformations
#' @examples
#' \dontrun{
#' ########################################################################
#' ## Example 1: Steady-state trafo
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pi(f, "A")
#' 
#' p.outerValues <- c(k1 = 1, k2 = 0.1, A = 10, B = 1)
#' P.steadyState(p.outerValues)
#' 
#' ########################################################################
#' ## Example 2: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pi(f, "A")
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState %o% P.log)(p.outerValue)
#' }
#' @export
Pi <- function(trafo, parameters=NULL, compile = FALSE) {

  states <- names(trafo)
  nonstates <- getSymbols(trafo, exclude = states)
  dependent <- setdiff(states, parameters)
  
  trafo.alg <- funC0(trafo[dependent], compile = compile)
  ftrafo <- function(x, parms) {
    out <- trafo.alg(as.list(c(x, parms)))
    structure(as.numeric(out), names = colnames(out))
  }
  
  # the parameter transformation function to be returned
  p2p <- function(p, fixed=NULL, deriv = TRUE) {
    
    # Inherit from p
    dP <- attr(p, "deriv")
    
    # replace fixed parameters by values given in fixed
    if(!is.null(fixed)) {
      is.fixed <- which(names(p)%in%names(fixed)) 
      if(length(is.fixed)>0) p <- p[-is.fixed]
      p <- c(p, fixed)
    }
    
    # check for parameters which are not computed by multiroot
    emptypars <- names(p)[!names(p)%in%c(dependent, fixed)]
    
    # Compute steady state concentrations
    myroot <- rootSolve::multiroot(ftrafo, 
                                   start = p[dependent], 
                                   parms = p[setdiff(names(p), dependent)])
    
    # Output parameters
    out <- c(myroot$root, p[setdiff(names(p), names(myroot$root))])
    
    # Compute jacobian d(root)/dp
    dfdx <- rootSolve::gradient(ftrafo, 
                                x = myroot$root, 
                                parms = p[setdiff(names(p), names(myroot$root))])
    dfdp <- rootSolve::gradient(ftrafo, 
                                x = p[setdiff(names(p), names(myroot$root))],
                                parms = myroot$root)
    dxdp <- solve(dfdx, -dfdp)
    #print(dxdp)
        
       
    # Assemble total jacobian
    jacobian <- matrix(0, length(out), length(p))
    colnames(jacobian) <- names(p)
    rownames(jacobian) <- names(out)
    for(ep in emptypars) jacobian[ep, ep] <- 1
    jacobian[rownames(dxdp), colnames(dxdp)] <- dxdp 
    jacobian <- jacobian[,setdiff(names(p), names(fixed))]
    
    # Multiplication with deriv of p
    if(!is.null(dP)) jacobian <- jacobian%*%dP[colnames(jacobian),]
    
    
    if(deriv) attr(out, "deriv") <- jacobian
    return(out)
    
  }
  
  class(p2p) <- "par"
  
  attr(p2p, "equations") <- trafo
  attr(p2p, "parameters") <- parameters
  
  return(p2p)
  
}


#' Concatenation of parameter transformations
#' 
#' @param p1 Return value of \link{P} or \link{Pi}
#' @param p2 Return value of \link{P} or \link{Pi}
#' @return A function \code{p2p(p, fixed = NULL, deriv = TRUE)}, the concatenation of \code{p1} and 
#' \code{p2}.
#' @aliases concatenation
#' @examples
#' \dontrun{
#' #' ########################################################################
#' ## Example: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pi(f, "A")
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState %o% P.log)(p.outerValue)
#' }
#' @export
"%o%" <- function(p1, p2) function(p, fixed=NULL, deriv = TRUE) p1(p2(p, fixed = fixed, deriv = deriv), deriv = deriv)

