
plaqr <- function(formula, nonlinVars=NULL, tau=.5, data=NULL, subset, basis=3,  
            weights, na.action, method = "br", model = TRUE, contrasts = NULL, ...)
{
    if(length(tau)>1) stop("tau must be a number (not a vector) strictly between 0 and 1.")
    if(class(nonlinVars)!="formula" & !is.null(nonlinVars)){
      stop("nonlinVars must be of class \"formula\" or NULL. NULL by default.\n")
    }
    plaqrcall <- match.call()
    rqcall <- plaqrcall
    nonlinvars <- NULL
    linvars <- attr(terms(formula, data=data), "term.labels")

    if(is.null(nonlinVars)){
      int <- ifelse(attr(terms(formula, data=data),"intercept")==1, "1", "0")
      rqcall$formula <- update(formula, paste(c("~",linvars,int),
                             collapse="+"))
    } else {
      nonlinvars <- attr(terms(nonlinVars, data=data), "term.labels")
      nonlinvars <- nonlinvars[!(nonlinvars %in% all.vars(formula)[1])]
      linvars <- linvars[!(linvars %in% nonlinvars)]
      nonlinvarsbs <- apply(matrix(nonlinvars), 1,
                         function(x) paste("bs(",x,",degree=",basis,")",sep=""))
      rqcall$formula <- update(formula, paste(c("~","1",linvars,nonlinvarsbs),
                             collapse="+"))
    }

    rqcall[[1]] <- as.name("rq")
    rqcall$nonlinVars <- rqcall$basis <- NULL
    model <- eval.parent(rqcall)
    class(model) <- c("plaqr", "rq")
    model$call <- plaqrcall
    model$linear <- linvars
    model$nonlinear <- nonlinvars
    if(is.null(nonlinVars)){
      model$z <- data.frame()
    } else{
      model$z <- model.frame(nonlinVars, data=data)
    }
    model$basis <- basis
    if(length(model$linear)==0) model$linear <- vector("character", 0)
    if(length(model$nonlinear)==0) model$nonlinear <- vector("character", 0)

    return(model)
}
