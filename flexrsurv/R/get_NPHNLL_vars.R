all_specials_vars.terms <- function(termsobj,
                                  specials="NPHNLL",
                                  functions = FALSE,
                                  max.names = 1L,
                                  unique = TRUE,
                                  order=c("specials", "formula")){
# extracts the names of the variables in specials effects from a terms object
# output : vector of character names  
  order <- match.arg(order)
  indxvars <- unlist(attr(termsobj, "specials")[specials])
  oindxvars <- order(indxvars)
  if(order=="specials"){
    allvars <- lapply(attr(termsobj,"term.labels")[indxvars-1]
                                        #       , function(x){all.vars(reformulate(x), max.names=1)})
                      , function(x){all.vars(reformulate(x),
                                             functions = functions,
                                             max.names = max.names,
                                             unique = unique)})
  }
  else {
    allvars <- lapply(attr(termsobj,"term.labels")[indxvars[oindxvars]-1]
#       , function(x){all.vars(reformulate(x), max.names=1)})
                      , function(x){all.vars(reformulate(x),
                                             functions = functions,
                                             max.names = max.names,
                                             unique = unique)})
  }    
  if(unique){
    unique(unlist(allvars))
  }
  else{
   unlist(allvars)
  }
}




get_all_specials_vars <- function(formula,
                                  data = NULL,
                                  specials="NPHNLL",
                                  max.names=1L,
                                  unique = TRUE,
                                  order=c("formula", "specials"),
                                  ...)
{
  # derived from stats::get_all_vars (2.15.1)
  # get all vars in NPHNLL terms
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
    varnames <- all_specials_vars.terms(formula, specials=specials, max.names=max.names,
                                unique = unique, order=order )
    inp <- parse(text=paste("list(", paste(varnames, collapse=","), ")"))
    variables <- eval(inp, data, env)
    if(is.null(rownames) && (resp <- attr(formula, "response")) > 0) {
        ## see if we can get rownames from the response
        lhs <- variables[[resp]]
        rownames <- if(is.matrix(lhs)) rownames(lhs) else names(lhs)
    }
    extras <- substitute(list(...))
    extranames <- names(extras[-1L])
    extras <- eval(extras, data, env)
    x <- as.data.frame(c(variables, extras), optional=TRUE)
    names(x) <- c(varnames, extranames)
    if (!is.null(rownames))
	attr(x, "row.names") <- rownames # might be short form
    x
}







