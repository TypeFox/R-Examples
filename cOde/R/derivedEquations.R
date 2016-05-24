
#' Compute sensitivity equations of a function symbolically
#' 
#' @param f named vector of type character, the functions
#' @param states Character vector. Sensitivities are computed with respect to initial
#' values of these states
#' @param parameters Character vector. Sensitivities are computed with respect to initial
#' values of these parameters
#' @param inputs Character vector. Input functions or forcings. They are excluded from
#' the computation of sensitivities.
#' @param reduce Logical. Attempts to determine vanishing sensitivities, removes their
#' equations and replaces their right-hand side occurences by 0.
#' @details The sensitivity equations are ODEs that are derived from the original ODE f.
#' They describe the sensitivity of the solution curve with respect to parameters like 
#' initial values and other parameters contained in f. These equtions are also useful
#' for parameter estimation by the maximum-likelihood method. For consistency with the
#' time-continuous setting provided by \link{adjointSymb}, the returned equations contain
#' attributes for the chisquare functional and its gradient.
#' @return Named vector of type character with the sensitivity equations. Furthermore,
#' attributes "chi" (the integrand of the chisquare functional), "grad" (the integrand
#' of the gradient of the chisquare functional), "forcings" (Character vector of the 
#' additional forcings being necessare to compute \code{chi} and \code{grad}) and "yini" (
#' The initial values of the sensitivity equations) are returned.
#' @example inst/examples/example2.R
#' @example inst/examples/example3.R
#' @export
sensitivitiesSymb <- function(f, states = names(f), parameters = NULL, inputs = NULL, reduce = FALSE) {
  
  variables <- names(f)
  states <- states[!states%in%inputs]
  
  if(is.null(parameters)) {
    pars <- getSymbols(f, exclude=c(variables, inputs, "time"))
  } else {
    pars <- parameters[!parameters%in%inputs]
  }
  
  Dyf <- jacobianSymb(f, variables)
  Dpf <- jacobianSymb(f, pars)
  
  df <- length(f)
  dv <- length(variables)
  ds <- length(states)
  dp <- length(pars)
  
  # generate sensitivity variable names and names with zero entries
  mygridY0 <- expand.grid.alt(variables, states)
  mygridP <- expand.grid.alt(variables, pars)
  sensParVariablesY0 <- apply(mygridY0, 1, paste, collapse = ".")
  sensParVariablesP <- apply(mygridP, 1, paste, collapse = ".")
  
  # Write sensitivity equations in matrix form
  Dy0y <- matrix(sensParVariablesY0, ncol = ds, nrow = dv)
  Dpy <- matrix(sensParVariablesP, ncol = dp, nrow = dv)
  
  
  gl <- c(as.vector(prodSymb(matrix(Dyf, ncol = dv), Dy0y)), 
          as.vector(sumSymb(prodSymb(matrix(Dyf, ncol = dv), Dpy), matrix(Dpf, nrow = dv))))
 
  newfun <- gl
  newvariables.grid <- expand.grid.alt(variables, c(states, pars))
  newvariables <- apply(newvariables.grid, 1, paste, collapse=".")
  names(newfun) <- newvariables
  
  # Reduce the sensitivities
  vanishing <- c(sensParVariablesY0, sensParVariablesP[Dpf == "0"])
  if(reduce) {
    newfun <- reduceSensitivities(newfun, vanishing)
    is.zero.sens <- names(newfun) %in% attr(newfun,"is.zero")
  } else {
    is.zero.sens <- rep(FALSE, length(newfun))
  }
  newfun <- newfun[!is.zero.sens]
  
  
    
  # Append initial values
  initials <- rep(0, length(newfun))
  names(initials) <- newvariables[!is.zero.sens]
  ones <- which(apply(newvariables.grid, 1, function(row) row[1] == row[2]))
  initials[newvariables[ones]] <- 1
  
  # Compute wrss
  pars <- c(pars, states)
  
  statesD <- paste0(states, "D")
  weightsD <- paste0("weight", statesD)
  
  res <- paste0(weightsD,"*(", states, "-", statesD, ")")
  sqres <- paste0(res, "^2")
  chi <- c(chi = paste0(sqres, collapse =" + "))
    
  sensitivities <- lapply(pars, function(p) paste0(states, ".", p))
  names(sensitivities) <- pars
    
  grad <- sapply(pars, function(p) paste0(paste("2*", res, "*", sensitivities[[p]]), collapse=" + "))
  names(grad) <- paste("chi", pars, sep=".")
  grad <- replaceSymbols(newvariables[is.zero.sens], "0", grad)
  
  attr(newfun, "chi") <- chi
  attr(newfun, "grad") <- grad
  attr(newfun, "outputs") <- structure(rep(0, length(which(is.zero.sens))), names = newvariables[is.zero.sens])
  attr(newfun, "forcings") <- c(statesD, weightsD)
  attr(newfun, "yini") <- initials
    
  
  return(newfun)
  
  
}




#' Compute adjoint equations of a function symbolically
#' 
#' @param f Named vector of type character, the functions
#' @param states Character vector of the ODE states for which observations are available
#' @param inputs Character vector of the "variable" input states, i.e. time-dependent parameters
#' (in contrast to the forcings).
#' @param parameters Character vector of the parameters
#' @details The adjoint equations are computed with respect to the functional 
#' \deqn{(x, u)\mapsto \int_0^T \|x(t)-x^D(t)\|^2 + \|u(t) - u^D(t)\|^2 dt,}{(x, u) -> int( ||x(t) - xD(t)||^2 + ||u(t) - uD(t)||^2, dt),} 
#' where x are the states being constrained
#' by the ODE, u are the inputs and xD and uD indicate the trajectories to be best
#' possibly approached. When the ODE is linear with respect to u, the attribute \code{inputs}
#' of the returned equations can be used to replace all occurences of u by the corresponding
#' character in the attribute. This guarantees that the input course is optimal with
#' respect to the above function.
#' @return Named vector of type character with the adjoint equations. The vector has attributes
#' "chi" (integrand of the chisquare functional), "grad" (integrand of the gradient of the chisquare functional),
#' "forcings" (character vector of the forcings necessary for integration of the adjoint equations) and
#' "inputs" (the input expressed as a function of the adjoint variables).
#' @example inst/examples/example5.R
#' @export
adjointSymb <- function(f, states=names(f), parameters = NULL, inputs=NULL) {
  
  n <- length(f)
  
  adjNames <- paste0("adj", names(f))
  
  ## Compute adjoint sensitivities  
  jac <- matrix(jacobianSymb(f), n)
  
  negadj <- as.vector(prodSymb(t(jac), matrix(adjNames, ncol=1)))
  adj <- paste0("-(", negadj, ")")
  names(adj) <- adjNames
  
  
  negres <- paste0("(", states, "D - ", states, ")") 
  wres <- paste(negres, paste0("weight", states, "D"), sep="*") 
  
  
  adj[paste0("adj", states)] <- paste(adj[paste0("adj", states)], wres, sep="+")
  
  ## Input equations
  u <- NULL
  if(!is.null(inputs)) {
    
    jac <- matrix(jacobianSymb(f, inputs), n )
    u <- as.vector(
      sumSymb(paste0(inputs, "D"), 
              matrix(paste0("-(",  
                            as.vector(
                              prodSymb(t(jac), matrix(adjNames, ncol=1))),
                            ")*eps/weight", inputs, "D"),
                     ncol=1)
      )
    )
    u <-  paste0("(", u, ")")
    names(u) <- inputs
    
  }
  
  ## Forcings required by the BVP solver
  forcings <- paste0(c(states, inputs), "D")
  forcings <- c(forcings, paste0("weight", forcings))
  
  
  
  ## Time-continous log-likelihood
  res <- paste0("(", c(states, inputs), "-", c(states, inputs), "D)")
  sres <- paste0(res, "^2")
  wsres <- paste(paste0("weight", c(states, inputs), "D"),sres,  sep="*") 
  
  chi <- paste(wsres, collapse = " + ")
  names(chi) <- "chi"
  attr(adj, "chi") <- chi
  
  ## Adjoint "gradient wrt parameters" equations
  if(is.null(parameters)) {
    symbols <- getSymbols(f)
    parameters <- symbols[!(symbols%in%inputs) & !(symbols%in%names(f))]
  }
  if(length(parameters)>0) {
    jac <- matrix(jacobianSymb(f, parameters), n)
    grad <- as.vector(prodSymb(matrix(adjNames, nrow=1), jac))
    gradP <- paste0("2*(", grad, ")")
    names(gradP) <- paste("chi", parameters, sep=".")
    attr(adj, "grad") <- gradP  
  }
  
  
  
  attr(adj, "forcings") <- forcings
  attr(adj, "inputs") <- u
  
  
  
  
  return(adj)
  
}


#'@title Forcings data.frame
#'@name forcData
#'@docType data
#' 
NULL

#'@title Time-course data of O, O2 and O3
#'@name oxygenData
#'@docType data 
NULL
