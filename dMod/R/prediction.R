#' Model evaluation. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function \code{x(times, pars, forcings, events, deriv = TRUE)} returning ODE output and sensitivities.
#' @param func return value from \code{funC(f)} where \code{f} defines the ODE.
#' @param extended return value from \code{funC(c(f, sensitivitiesSymb(f)))}.
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link[deSolve]{events}.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @return A model prediction function \code{x(times, pars, forcings, events, deriv = TRUE)} representing 
#' the model evaluation. The result of
#' \code{x(times, pars, forcings, events, deriv = TRUE)} contains
#' attributes "sensitivities" and "deriv" with the sensitivities if \code{deriv=TRUE}. 
#' If \code{deriv=FALSE}, sensitivities are not computed (saving time).
#' If \code{pars} is
#' the result of \code{p(pouter)} (see \link{P}), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv", 
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide. 
#' @export
Xs <- function(func, extended, forcings=NULL, events=NULL, optionsOde=list(method="lsoda"), optionsSens=list(method="lsodes")) {
  
  myforcings <- forcings
  myevents <- events
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables")%in%variables]
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) v[2]))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  
  
  ## Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar

  # Names for deriv output
  sensGrid <- expand.grid(variables, c(svariables, sparameters), stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  P2X <- function(times, pars, forcings = myforcings, events = myevents, deriv=TRUE){
    
    myforcings <- forcings
    myevents <- events
    
    yini <- pars[variables]
    mypars <- pars[parameters]
    
    
    if(!deriv) {
    
      # Evaluate model without sensitivities
      loadDLL(func)
      if(!is.null(myforcings)) forc <- setForcings(func, myforcings) else forc <- NULL
      out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=mypars, forcings=forc, events = list(data = events)), optionsOde))
      #out <- cbind(out, out.inputs)
      
      
    } else {
      
      # Evaluate extended model
      loadDLL(extended)
      if(!is.null(myforcings)) forc <- setForcings(extended, myforcings) else forc <- NULL
      outSens <- do.call(odeC, c(list(y=c(yini, yiniSens), times=times, func=extended, parms=mypars, forcings=forc, events = list(data = events)), optionsSens))
      #out <- cbind(outSens[,c("time", variables)], out.inputs)
      out <- outSens[,c("time", c(variables, forcnames))]
      attr(out, "sensitivities") <- outSens[,!colnames(outSens)%in%c(variables, forcnames)]
      
      
      # Apply parameter transformation to the derivatives
      sensLong <- matrix(outSens[,sensNames], nrow=dim(outSens)[1]*length(variables))
      dP <- attr(pars, "deriv")
      if(!is.null(dP)) {
        sensLong <- sensLong%*%(dP[c(svariables, sparameters),])
        sensGrid <- expand.grid(variables, colnames(dP), stringsAsFactors=FALSE)
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      }
      outSens <- cbind(outSens[,1], matrix(sensLong, nrow=dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      
      attr(out, "deriv") <- outSens
      attr(out, "parameters") <- unique(sensGrid[,2])
    }
    
    return(out)
    
  }
  
  return(P2X)
  
  
}



#' Observation functions. 
#' @description Creates a function \code{y(out, pars)} that evaluates an observation function
#' and its derivatives based on the output of a model function \code{x(times, pars)}, see \link{Xf} and \link{Xs}.
#' @param g Named character vector defining the observation function
#' @param f Named character, the underlying ODE
#' @param compile Logical, compile the function (see \link{funC0})
#' @return a function \code{y(out, pars, attach=FALSE)} representing the evaluation of the observation function. 
#' If \code{out} has the attribute  "sensitivities", the result of
#' \code{y(out, pars)}, will have an attributed "deriv" which reflec the sensitivities of 
#' the observation with respect to the parameters.
#' If \code{pars} is the result of a parameter transformation \code{p(pars)} (see \link{P}), 
#' the Jacobian 
#' of the parameter transformation and the sensitivities of the observation function
#' are multiplied according to the chain rule for differentiation.
#' If \code{attach = TRUE}, the original argument \code{out} will be attached to the evaluated observations.
#' @export
Y <- function(g, f, compile = FALSE) {
  
  warnings <- FALSE
  
  # Get potential paramters from g, forcings are treated as parameters because
  # sensitivities dx/dp with respect to forcings are zero.
  states <- names(f)
  parameters <- getSymbols(c(g, f), exclude = c(states, "time"))
    
  # Observables defined by g
  observables <- names(g)
  gEval <- funC0(g, compile = compile)
  
  # Character matrices of derivatives
  dxdp <- dgdx <- dgdp <- NULL
  
  if(length(states) > 0 & length(parameters) > 0) {
    dxdp <- apply(expand.grid.alt(states, c(states, parameters)), 1, paste, collapse = ".")
    dxdp <- matrix(dxdp, nrow = length(states))
  }
  if(length(states) > 0)
    dgdx <- matrix(jacobianSymb(g, states), nrow=length(g))
  if(length(parameters) > 0) {
    dgdp <- cbind(
      matrix("0", nrow=length(g), ncol=length(states)), 
      matrix(jacobianSymb(g, parameters), nrow=length(g))
    )
  }
  
  # Sensitivities of the observables
  derivs <- as.vector(sumSymb(prodSymb(dgdx, dxdp), dgdp))
  if(length(derivs) == 0) stop("Nor states or parameters involved")
  names(derivs) <- apply(expand.grid.alt(observables, c(states, parameters)), 1, paste, collapse = ".")
  derivsEval <- funC0(derivs, compile = compile)
    
  # Vector with zeros for possibly missing derivatives
  zeros <- rep(0, length(dxdp))
  names(zeros) <- dxdp
  
  
  X2Y <- function(out, pars, attach=FALSE) {
    
    # Prepare list for with()
    nOut <- dim(out)[2]
    outlist <- lapply(1:nOut, function(i) out[,i]); names(outlist) <- colnames(out)
    
    dout <- attr(out, "sensitivities")
    if(!is.null(dout)) {
      nDeriv <- dim(dout)[2]
      derivlist <- lapply(1:nDeriv, function(i) dout[,i]); names(derivlist) <- colnames(dout)  
    } else {
      derivlist <- NULL
    }
    
      
    x <- c(outlist, derivlist, as.list(pars), as.list(zeros))
    
    values <- gEval(x)
    if(!is.null(dout)) dvalues <- derivsEval(x)
  
    # Parameter transformation
    dP <- attr(pars, "deriv")
    if(!is.null(dP) & !is.null(dout)) {
      
      parameters.all <- c(states, parameters)
      parameters.missing <- parameters.all[!parameters.all%in%rownames(dP)]
      
      if(length(parameters.missing) > 0 & warnings)
        warning("Parameters ", paste(parameters.missing, collapse = ", ", "are missing in the Jacobian of the parameter transformation. Zeros are introduced."))
      
      dP.missing <- matrix(0, nrow = length(parameters.missing), ncol=dim(dP)[2], 
                           dimnames=list(parameters.missing, colnames(dP)))
      dP <- rbind(dP, dP.missing)
      
      # Multiplication with tangent map
      sensLong <- matrix(dvalues, nrow=dim(out)[1]*length(observables))
      sensLong <- sensLong%*%(dP[parameters.all,])
      dvalues <- matrix(sensLong, nrow=dim(out)[1])
      
      # Naming
      sensGrid <- expand.grid.alt(observables, colnames(dP))
      sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      colnames(dvalues) <- sensNames
      
    }
  
  
  
    # Format output
    values <- cbind(time = out[,"time"], values)
    if(attach) 
      values <- cbind(values, out[,-1])
    
    
    if(!is.null(dout) && !attach) {
      attr(values, "deriv") <- cbind(time = out[,"time"], dvalues)
      if(is.null(dP)) attr(values, "parameters") <- names(pars) else attr(values, "parameters") <- colnames(dP)
    }
    if(!is.null(dout) && attach) {
      attr(values, "deriv") <- cbind(time = out[,"time"], dvalues, attr(out, "deriv")[,-1])
      if(is.null(dP)) attr(values, "parameters") <- names(pars) else attr(values, "parameters") <- colnames(dP)
    }
    return(values)        
  

  }


  attr(X2Y, "equations") <- g
  return(X2Y)


}

#' Model evaluation without sensitivities. 
#' @description Interface to get an ODE 
#' into a model function \code{x(times, pars, forcings, events)} returning ODE output.
#' It is a reduced version of \link{Xs}, missing the sensitivities. 
#' @param func return value from \code{funC(f)} where \code{f} defines the ODE. 
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link[deSolve]{events}.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @details Can be used to integrate additional quantities, e.g. fluxes, by adding them to \code{f}. All quantities that are not initialised by pars 
#' in \code{x(times, pars, forcings, events)} are initialized at 0.
#' @export
Xf <- function(func, forcings=NULL, events=NULL, optionsOde=list(method="lsoda")) {
  
  myforcings <- forcings
  myevents <- events
  
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0,length(variables))
  names(yini) <- variables
  
  P2X <- function(times, P, changedForcings = NULL, events = myevents){
    
    if(!is.null(changedForcings)) myforcings <- changedForcings
    yini[names(P[names(P) %in% variables])] <- P[names(P) %in% variables]
    pars <- P[parameters]
    #alltimes <- unique(sort(c(times, forctimes)))
    
    loadDLL(func)
    if(!is.null(myforcings)) forc <- setForcings(func, myforcings) else forc <- NULL
    out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=pars, forcings=forc,events = list(data = events)), optionsOde))
    #out <- cbind(out, out.inputs)      
    
    return(out)
    
  }
  
  return(P2X)
  
}
