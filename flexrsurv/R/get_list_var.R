get_all_Xvars.formula <- function(formula, unique=TRUE){
# output : the names of the X variables in the formula

  allvars <- all.vars(update(formula, 1 ~.))
    if(unique){
      unique(allvars)
    }
    else{
      allvars
    }
}
  
all_LIN_vars<- function(objterm){
  indxspc <- unlist(attr(objterm, "specials"))
  return(as.character(attr(objterm,"variables"))[-c(1,2, indxspc+1)])
}


all_specials_vars <- function(objterm, specials=c("NLL", "NPH", "NPHNLL"),
                              unique = TRUE,
                              order=c("formula", "specials"),
                            ...){
# extracts the names of the variables in specials effects from a terms object
# output : vector of character names  

  order <- match.arg(order)

  indxvar <- unlist(attr(objterm, "specials")[specials])
  nvars <- length(indxvar)

#  print(c(nvars, -1, indxvar))
  if(nvars==0){
    # no "specials" vars 
    return(character(0))
  }
  else{
    if(order=="specials"){
      oindxvar <- 1:nvars
    }
    else {
      oindxvar <- order(indxvar)
    }  
    var_list <- NULL
    for( i in indxvar[oindxvar]){

      fun <- mget(as.character(attr(objterm,"variables")[[i+1]][[1]]),
                  mode = "function",
                  envir = parent.frame(), inherits=TRUE,
                  ifnotfound=list(NULL))[[1]]
      if( !is.null(fun) ){
        thecall <- match.call(fun, attr(objterm,"variables")[[i+1]])
# the variable name is the second argument of the spetial function
        var_list <- c( var_list, thecall[[2]])
      }
    }
    if(unique){
      as.character(unique(unlist( var_list)))
    }
    else{
      as.character(unlist(var_list))
    }
  }

}

all_NPHNLL_timevar <- function(objterm, specials=c("NPH", "NPHNLL"),
                              unique = FALSE,
                              order=c("formula", "specials"),
                               ...){

  order <- match.arg(order)

  indxvar <- unlist(attr(objterm, "specials")[specials])

  
  nvars <- length(indxvar)
  if(nvars==0){
    # no "specials" vars 
    return(character(0))
  }
  else{
    if(order=="specials"){
      oindxvar <- 1:nvars
    }
    else {
      oindxvar <- order(indxvar)
    }  
    var_list <- NULL
    for( i in indxvar[oindxvar]){
      fun <- mget(as.character(attr(objterm,"variables")[[i+1]][[1]]),
                  mode = "function",
                  envir = parent.frame(), inherits=TRUE,
                  ifnotfound=list(NULL))[[1]]
      if( !is.null(fun) ){
        thecall <- match.call(fun, attr(objterm,"variables")[[i+1]])
      # the time_variable name is the 3rd argument of the special function
        var_list <- c( var_list, thecall[[3]])
      }
    }

    if(unique){
      as.character(unique(unlist( var_list)))
    }
    else{
      as.character(unlist(var_list))
    }

  }

}





get_specials_vars <- function(formula,
                              data = NULL,
                              specials="NPHNLL",
                              all.vars.func=all_specials_vars, 
                              unique = TRUE,
                              order=c("formula", "specials"),
                              ...)
{
  # derived from stats::get_all_vars (2.15.1)
  # get vars selected by the function all.vars.func

  
    if(missing(formula)) {
	if(!missing(data) && inherits(data, "data.frame") &&
	   length(attr(data, "terms")) )
	    return(data)
	formula <- as.formula(data)
    }
    else if(missing(data) && inherits(formula, "data.frame")) {
	if(length(attr(formula, "terms")))
	    return(formula)
	data <- formula
	formula <- as.formula(data)
    }
    formula <- as.formula(formula)
    if(missing(data))
	data <- environment(formula)
    else if (!is.data.frame(data) && !is.environment(data)
             && !is.null(attr(data, "class")))
        data <- as.data.frame(data)
    else if (is.array(data))
        stop("'data' must be a data.frame, not a matrix or an array")
    if(!inherits(formula, "terms"))
	formula <- terms(formula, data = data)
    env <- environment(formula)
    rownames <- .row_names_info(data, 0L) #attr(data, "row.names")
# cat("    # colnames data\n")
#    print(data[[1]][1:20])
#cat("    # colnames data\n")
    # choose the variables
#cat("\n    # choose the variables\n")
    varnames <- all.vars.func(formula, specials=specials,  
                              order=order, unique = unique, ... )
#cat("    # choose the variables\n")
    inp <- parse(text=paste("list(", paste(varnames, collapse=","), ")"))
#cat(as.character(inp), "\n    variables <- eval(inp, data, env)\n")
    variables <- eval(inp, data, env)
#cat("    variables <- eval(inp, data, env)\n")
    if(is.null(rownames) && (resp <- attr(formula, "response")) > 0) {
        ## see if we can get rownames from the response
        lhs <- variables[[resp]]
        rownames <- if(is.matrix(lhs)) rownames(lhs) else names(lhs)
    }
    extras <- substitute(list(...))
    extranames <- names(extras[-1L])
 #   print(extras)
 #   print(extranames)
#cat("    extras <- eval(extras, data, env)\n")
    extras <- eval(extras, data, env)
#cat("    extras <- eval(extras, data, env)\n")
    x <- as.data.frame(c(variables, extras), optional=TRUE)
    names(x) <- c(varnames, extranames)
    if (!is.null(rownames))
	attr(x, "row.names") <- rownames # might be short form
    x
}


