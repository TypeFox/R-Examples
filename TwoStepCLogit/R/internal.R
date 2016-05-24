# Internal functions in the package TwoStepCLogit.
# These are not documented and they are hidden in the namespace.

shared <- function(formula, data, data.name){

  call <- match.call()
  
  # Validation de l'argument data
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (grepl("(", data.name, fixed=TRUE)) stop("the 'data' argument must be an object, not a function call", call. = FALSE)
  
  # Validation de l'argument 'formula'
  if (class(formula)!="formula") stop("the 'formula' argument must be a formula", call. = FALSE)
  ## covariables
  formula_terms <- terms(formula, special=c("strata", "cluster"), data = data)
  info.cluster <- untangle.specials(formula_terms, 'cluster')
  if( length(info.cluster$vars) == 0) stop("the formula must include a cluster term", call. = FALSE)
  info.strata <- untangle.specials(formula_terms, 'strata')
  if( length(info.strata$vars) == 0) stop("the formula must include a strata term", call. = FALSE)
  ## variable reponse
  appel.update_nocs <- paste("update.formula(old=formula, new= ~ . -", info.cluster$vars, "-", info.strata$vars, ")") 
  formula_nocs <- eval(parse(text=appel.update_nocs))
  mf <- model.frame(formula_nocs, data=data, na.action=NULL, drop.unused.levels=TRUE)
  y <- model.response(mf)
  if (is.null(y)) stop("a response variable must be given", call. = FALSE)
  if (!all(y %in% 0:1)) stop("the response variable can only take the values 0 and 1", call. = FALSE)
  mm <- model.matrix(formula_nocs, data=mf)
  
  # Pour retirer l'information sur les clusters
  var.cluster <- eval(parse(text=info.cluster$vars), envir=data)
  if (grepl("+", info.cluster$vars, fixed=TRUE) || length(dim(var.cluster)) > 1) 
    stop("in 'formula', the cluster identifier must be a single variable", call. = FALSE)
  
  # Sortie des resultats
  out <- list(data=data, info.cluster=info.cluster, info.strata=info.strata, y=y, mm=mm, var.cluster=var.cluster)
  return(out)
}


tryCatch.W.E <- function(expr)
{ # Fonction pour stocker les erreurs et les warnings
  # tiree de demo(error.catching), legerement modifiee
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W, w$message) ## ma modif ici pour stocker tous les warnings plutot que seulement le dernier
    invokeRestart("muffleWarning")
  }
  e.handler <- function(e){ # error handler
    class(e) <- "erreur"
    return(e)
  }
  list(value = withCallingHandlers(tryCatch(expr, error = e.handler),
          warning = w.handler),
      warnings = W)
}

