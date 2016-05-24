# nocov start
# ************************************************************************************
# TO DO: 
# REIMPORT THIS FILE BACK INTO SIMCAUSAL
# ************************************************************************************
  # x) Extend the checking for non_TD pars to TD parents in find_FormVars() 
  # x) For non_TD var outside of the DAG also check that length(non_TD) < 2
  # x) => Want to allow vectors in user.env to be referenced as uservec[t]
  # x) => This will allow avoiding declaration of node attributes as nodes, will save a ton of memory
# ************************************************************************************

opts <- new.env(parent = emptyenv())
opts$NoChangeFunCalls <- TRUE # Flag, if TRUE will not modify any unknown node formula functions while parsing
opts$vecfun <- NULL           # character vector of user-defined vectorized function names
opts$debug <- FALSE            # debug mode, when TRUE print all calls to dprint()

dprint <- function(...) if (opts$debug) print(...) # debug-only version of print
debug_set <- function() { # Set to Debug Mode
  mode <- TRUE
  old <- opts$debug
  opts$debug <- mode
  invisible(old)
}
debug_off <- function() { # Turn Off Debug Mode
  mode <- FALSE
  old <- opts$debug
  opts$debug <- mode
  invisible(old)
}
vecfun.print <- function() {
  new <- opts$vecfun
  if (length(new)>1) new <- paste0(new, collapse=",")
  print("current list of user-defined vectorized functions: "%+%new)
  invisible(opts$vecfun)
}
vecfun.all.print <- function() {
  new <- opts$vecfun
  if (length(new)>1) new <- paste0(new, collapse=",")
  print("build-in vectorized functions:"); print(c(vector_ops_fcns, vector_math_fcns))
  print("user-defined vectorized functions: "%+%new)
  invisible(c(vector_ops_fcns,vector_math_fcns,new))
}
vecfun.add <- function(vecfun_names) { # Add vectorized function to global custom vectorized function list and return the old version
  old <- opts$vecfun
  opts$vecfun <- unique(c(opts$vecfun, vecfun_names))
  new <- opts$vecfun
  if (length(new)>1) new <- paste0(new, collapse=",")
  print("current list of user-defined vectorized functions: "%+%new)
  invisible(old)
}
vecfun.remove <- function(vecfun_names) { # Remove vectorized functions to global custom vectorized function list and return the old version
  old <- opts$vecfun
  idx_remove <- old%in%vecfun_names

  if (sum(idx_remove) < length(vecfun_names)) {
    fun_notfound <- vecfun_names[!(vecfun_names%in%old)]
    if (length(fun_notfound)>1) fun_notfound <- paste0(fun_notfound, collapse=",")
    warning("some of the function names in 'vecfun_names' were not found and cannot be removed: "%+%fun_notfound)
  }
  if (sum(idx_remove)>0) {
    opts$vecfun <- opts$vecfun[-(which(idx_remove))]
  }
  new <- opts$vecfun
  if (length(new)>1) new <- paste0(new, collapse=",")
  print("current list of user-defined vectorized functions: "%+%new)
  invisible(old)
}
vecfun.reset <- function() {
  old <- opts$vecfun
  opts$vecfun <- NULL
  invisible(old)
}
vecfun.get <- function() opts$vecfun
get_opts <- function() opts$debug # Return Current Debug Mode Setting

# ------------------------------------------------------------------------------------------------------
# FUNCTIONS FOR PARSING AND EVALUATING NODE FORMULAS (PARAMETERS)
# ------------------------------------------------------------------------------------------------------

# FUNCTION NAMES THAT PRODUCE A VECTOR WHEN APPLIED TO A VECTOR (WILL NOT BE REPLACED BY apply(df,1,func)):
vector_fcns <- c("cbind_mod","vecapply","apply","rowSums","rowMeans", "(", "[", "[[", "{", ":", "rep", "length", "if")
# vectorized operators:
vector_ops_fcns <- c("ifelse", "+", "-", "*","^", "/", "==", "!=", "!", "<", ">", "<=", ">=", "|", "&")
# vectorized math funcs
vector_math_fcns <- c("I","abs","sign","sqrt","round","signif","floor","ceil","ceiling","trunc",
                      "sin","tan","cos","acos","asin","atan","cosh","sinh","tanh",
                      "log","log10","log1p","exp","expm1","plogis",
                      "beta","lbeta","gamma","lgamma","psigamma","digamma","trigamma",
                      "choose","lchoose","factorial","lfactorial")

# a) find TD var calls;
# b) find baseline var calls;
# c) parse the tree at most 10 times and evaluate all atomic expressions
# d) modify calls to summary (non-vectorized) function to apply(DF, 1, func_name), adding cbind to calls with more than 1 arg
nodeform_parsers = function(node_form_call, data.env, user.env)  {
  # combine all default vectorized funs + the user-specified vectorized function in global :
  vector_fcns_all <- c(vector_fcns, vector_ops_fcns, vector_math_fcns, vecfun.get())
  curr.dfvarnms <- data.env[["ANCHOR_ALLVARNMS_VECTOR_0"]]
  # 
  # (not USED) SUMMARY FCNS (non-vectorized): these will be always turned into apply(arg, 1, func)
  # summary_fcns <- c("c","all","any","sum","mean","prod","min","max","range")
  # (not USED) FOR FUTURE IMPLEMENTATION: FUNCTION NAMES THAT AREN'T ALLOWED IN FORMULA EXPRESSIONS:
  # banned_fcns <- c( "apply", "cbind", "&&", "||")

  # * recursively parse the call tree structure for a given expression, find call to '[' or a name, then output that name (TDVar name will be called as TDVar[])

  # ************************************************************************************
  # TO DO: 
  # Extend the same checks for non_TD var existance to TD vars => Want to allow vectors in user.env to be referenced as uservec[t]
  # When TDvar_t not in DAG, check that TD_var exists in user.env, check that its a vector and that length matches t length
  # Decide between method I & II for finding non_TD parents
  # Curently using method I for non_TD vars, plotting DAG will exclude 
  # ************************************************************************************

  find_FormVars <- function(x, vartype="TD") {
    if (is.name(x) & vartype=="non_TD") {
      dprint("x: " %+% as.character(x))

      # Method I: will find all variables referenced by node formula, including vars only defined in user.env and are not part of the DAG
      notis.fun <- eval(substitute(!is.function(try(get(as.character(x)),  silent = TRUE))), envir = data.env, enclos = user.env)
      dprint("is x not a fun? " %+% notis.fun)

      # Method II: will only identify vars that were defined in the DAG. Probably more stable.
      is.inDAG <- as.character(x) %in% curr.dfvarnms
      dprint("is x in DAG? " %+% is.inDAG);

      # CHECK FOR UNDECLARED VARS: Verify if x is even defined in user env if (!notis.fun & !is.inDAG)
      exists.x <- exists(as.character(x), where = user.env, inherits = FALSE) # exists.x <- exists(as.character(x), where = user.env, inherits = TRUE)
      dprint("does x exist in user.env? " %+% exists.x)

      # CHECK THAT ITS NOT A SPECIAL (RESERVED) VAR (nF)
      specialVar <- "nF"
      special <- as.character(x)%in%specialVar
      dprint("x is special? " %+% special)

      if (notis.fun && !is.inDAG && !exists.x && !special) stop("Undefined variable: " %+% as.character(x), call. = FALSE)

      # ****************************
      # *) For non_TD var outside of the DAG also check that length(non_TD) < 2
      # ****************************

      # if (is.inDAG) varnames <- as.character(x) # METHOD I declares only vars that exist in the DAG as a parent
      if (notis.fun) varnames <- as.character(x) # METHOD II declares any nonfun var as a parent

    } else if (is.atomic(x) || is.name(x)) {
      character()
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`[`)) && is.name(x[[2]])) {
        if (vartype=="TD") {
          varnames <- as.character(x[[2]])
        } else if (vartype=="TD_t") {
          varnames <- as.character(as.character(x[[2]]) %+% "_" %+% eval(x[[3]]))
        } else if (vartype=="non_TD") {
          varnames <- character()
          x[[2]] <- NULL
        } else {
          stop("unrecognized variable type")
        }
      } else {
        varnames <- character()
      }
      if (length(x)>1 & vartype=="non_TD") x[[1]] <- NULL
      unique(c(varnames, unlist(lapply(x, find_FormVars, vartype))))
    } else if (is.pairlist(x)) {
      unique(unlist(lapply(x, find_FormVars, vartype)))
    } else {
      stop("Don't know how to handle type ", typeof(x), call. = FALSE)
    }
  }

  # * iteratively parse the call tree and evaluate all functions with atomic args until identical tree is returned
  eval_all_atomic <- function(expr) {
    eval_atomic <- function (x, where = parent.frame()) {
      if (is.atomic(x) || is.name(x)) {
        x	# Leave unchanged
      } else if (is.call(x)) {
        # reached '[', '[[' or 'c' functions, don't need to parse any deeper, return this subtree intact
        if (((identical(x[[1]], quote(`[`)) || identical(x[[1]], quote(`[[`))) && is.name(x[[2]])) || identical(x[[1]], quote(c))) {
          x # Leave unchanged
        } else {
          atomargs_test <- sapply(2:length(x), function(i) is.atomic(x[[i]]))
          # dprint("call: "%+%x[[1]]); for (i in (2:length(x))) dprint(x[[i]]); dprint("all atomic?: "%+%all(atomargs_test))
          if (!all(atomargs_test) | identical(x[[1]], quote(`{`))) { # 1) either one or more args are non-atomic, then continue parsing
            as.call(lapply(x, eval_atomic, where = where))
          } else {
            # or 2) all args are atomic - then evalute and return result
            # dprint("all args atomic, evaluated: "); dprint(eval(x))
            eval(x)
          }
        }
      } else if (is.pairlist(x)) {
        as.pairlist(lapply(x, eval_atomic, where = where))
      } else { # User supplied incorrect input
        stop("Don't know how to handle type ", typeof(x), call. = FALSE)
      }
    } # end of eval_atomic()

    dprint("expression before atomic pre-eval: "); dprint(expr)

    preveval_atom_call <- expr
    i <- 1
    samecall <- FALSE # flag for parsed tree being identical to the previous pre-parsed call tree

    # loop for max 10 iterations or when call tree is no longer changing:
    while ((i <= 10) & (!samecall)) {
      eval_atom_call <- eval_atomic(preveval_atom_call)
      samecall <- identical(eval_atom_call, preveval_atom_call)
      # dprint("-------------");
      # dprint(samecall); 
      # dprint(eval_atom_call); 
      # dprint("-------------")
      preveval_atom_call <- eval_atom_call
      i <- i + 1
    }
    # dprint("expression after atomic pre-eval: "); dprint(eval_atom_call)
    eval_atom_call
  }

  # * TO DO: MIGHT REMOVE THIS FUNCTION COMPLETELY, EITHER ASSUME ALL FUNs ARE ALREADY VECTORIZED OR PRE-TEST FUNS FOR VECTORIZATION AND RETURN ERROR IF FUN RETURNS A NON-VECTOR (SCALAR)
  # * TO ADD: if call tree starts with '{' need to process each argument as a separate call and return a list of calls instead
  # * modify the call tree with apply for non-vectorized (summary) functions, also adding cbind_mod() for calls with more than one arg
  modify_call <- function (x, where = parent.frame()) {
    if (is.atomic(x) & length(x)>1) {
      x <- parse(text = deparse(x, width.cutoff = 500))[[1]]
      modify_call(x, where = where)	# continue parsing recursively, turning result back into call
    }
    if (is.atomic(x) || is.name(x)) {
      if (is.atomic(x)) dprint("atomic: "%+%x)
      if (is.name(x)) dprint("name: "%+%x)
      x	# Leave unchanged
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`[`)) && is.name(x[[2]])) {	# reached '[' function, don't need to parse any deeper, return this subtree intact
        x
      } else if (identical(x[[1]], quote(`[[`)) && is.name(x[[2]])) { # reached '[[' function, same as above
        x
      } else if (as.character(x[[1]]) %in% vector_fcns_all)  {  # these functions are already vectorized (if given a vector, will return a vector)
        # dprint(paste0("vectorized func: ",as.character(x[[1]])))
        as.call(lapply(x, modify_call, where = where))	# continue parsing recursively, turning result back into call
      } else {	# non-vectorized fun needs to be wrapped in vecapply, with args combined as cbind(arg1,arg2,...) for more than one arg
        # dprint(paste0("non-vectorized func: ",as.character(x[[1]])))
        if (identical(x[[1]], quote(c))) {
          x[[1]] <- quote(cbind_mod)	# check if the function is 'c', in which case replace call with 'cbind_mod'
          as.call(lapply(x, modify_call, where = where))	# continue parsing recursively, turning result back into call
        } else if (identical(x[[1]], quote(sum))) {
          x[[1]] <- quote(rowSums)	# check if the function is 'sum', in which case replace call with 'colSums'
          as.call(lapply(x, modify_call, where = where))	# continue parsing recursively, turning result back into call
        } else if (identical(x[[1]], quote(mean))) {
          x[[1]] <- quote(rowMeans)	# check if the function is 'mean', in which case replace call with 'colMeans'
          as.call(lapply(x, modify_call, where = where))	# continue parsing recursively, turning result back into call
        } else if (identical(x[[1]], quote(structure))) {
          modify_call(as.call(x[[2]]), where = where)  # continue parsing recursively, turning result back into call          
          # as.call(lapply(x[[2]], modify_call, where = where))  # continue parsing recursively, turning result back into call          
        
        # OS 09/22/15: Adding new global option that prevents any modifications of node formulas, with a warning
        } else if (opts$NoChangeFunCalls) {
          message("Warning: function '" %+% deparse(x[[1]]) %+% "' will be called as is, even though it is not on the recognized vectorized functions list; use at your own risk!")
          as.call(lapply(x, modify_call, where = where))

        } else {
          nargs <- length(x)-1
          if (nargs > 1) { # IF NON-VECTORIZED func has more than one argument, combine all args into one with cbind_mod
          # dprint("several args: "%+%x[[1]])
            newargs <- "cbind_mod("%+%deparse(x[[2]], width.cutoff=500)
            for (i in (3:length(x))) newargs <- newargs%+%","%+%deparse(x[[i]], width.cutoff=500)
            for (i in (length(x)):3) x[[i]] <- NULL
            newargs <- newargs%+%")"
            newexp <- parse(text=newargs)[[1]]
            x[[2]] <- newexp
          }
          reparsed_chr <- "vecapply("%+%deparse(x[[2]], width.cutoff=500) %+% ", 1, " %+% deparse(x[[1]], width.cutoff=500) %+% ")"
          reparsed_call <- parse(text=reparsed_chr)[[1]]
          print(x)
          x[[1]] <- reparsed_call
          x[[2]] <- NULL
          modify_call(x[[1]], where = where)	# continue parsing recursively, turning result back into call
        }
      }
    } else if (is.pairlist(x)) {
      as.pairlist(lapply(x, modify_call, where = where))
    } else { # User supplied incorrect input
      stop("Don't know how to handle type ", typeof(x), call. = FALSE)
    }
  }

  # eval_atom_call <- node_form_call						      # don't evaluate any atomic expressions
  eval_atom_call <- eval_all_atomic(node_form_call)		# pre-evaluate all atomic expressions

  # Parses the formula and gets all the variable names referenced as [] or as.name==TRUE
  Vnames <- find_FormVars(eval_atom_call, vartype="non_TD")	# returns unique names of none TD vars that were called as VarName
  TD_vnames <- find_FormVars(eval_atom_call, vartype="TD")	# returns unique names TDVar that were called as TDVar[indx]
  TD_t_vnames <- find_FormVars(eval_atom_call, vartype="TD_t") # returns unique names TDVar_t that were called as TDVar[indx]

  dprint("Vnames: "); dprint(Vnames)
  dprint("TD_vnames: "); dprint(TD_vnames)
  dprint("TD_t_vnames: "); dprint(TD_t_vnames)

  modified_call <- modify_call(eval_atom_call) 			# parse current call and replace any non-vectorized function with apply call (adding cbind_mod if more than one arg)
  dprint("modified_call"); dprint(modified_call)

  return(list(Vnames = Vnames, TD_vnames = TD_vnames, TD_t_vnames = TD_t_vnames, modified_call = modified_call))
}

eval.nodeform.full <- function(expr_call, expr_str, self, data.env) {
  # traverse the node formula call, return TDvar & Var names (node parents) and modify subst_call to handle non-vectorized (summary functions):
  parse_res <- nodeform_parsers(node_form_call = expr_call, data.env = data.env, user.env = self$user.env)

  # set the local variables in the formula node to their character values:
  Vnames  <- parse_res$Vnames
  TD_vnames <- parse_res$TD_vnames
  TD_t_vnames <- parse_res$TD_t_vnames
  modified_call <- parse_res$modified_call # modified call that has any non-vectorized function replaced with apply call (with cbind for more than one arg)

  dprint("----------------------")
  dprint("node: "%+%self$cur.node$name)
  dprint("original expr as call:"); dprint(expr_call)
  dprint("final exprs:"); dprint(parse_res$modified_call)
  dprint("----------------------")

  df_varnms <- data.env$ANCHOR_ALLVARNMS_VECTOR_0

  # check this TD Var doesn't exist in parent environment if df has TD Var => TD Var is not time-dependent and reference TDVar[t] is incorrect - throw exception
  for (TDname in TD_vnames) {
    if (TDname%in%df_varnms) stop(paste0("reference ", TDname, "[...]", " at node ", self$cur.node$name, " is not allowed; node ", TDname," was defined as time-invariant"))
  }

  # NO LONGER NEED, this check is now performed in find_FormVars
  # for (Vname in Vnames) {
  #   if (!(Vname%in%df_varnms)) stop(paste0("formula at node ", self$cur.node$name, " cannot be evaluated; node ", Vname," is undefined"))
  # }

  data.env <- c(self$node_fun, data.env)

  if (is.call(modified_call) && identical(try(modified_call[[1]]), quote(`{`))) { # check for '{' as first function, if so, remove first func, turn call into a list of calls and do lapply on eval
    modified_call_nocurl <- modified_call[-1]
    evaled_expr <- try(lapply(X = modified_call_nocurl, FUN = eval, envir = data.env, enclos = self$user.env))
  } else {
    evaled_expr <- try(eval(modified_call, envir = data.env, enclos = self$user.env)) # eval'ing expr in the envir of data.df
  }

  if(inherits(evaled_expr, "try-error")) {
    stop("error while evaluating node "%+% self$cur.node$name %+%" formula: \n"%+%parse(text = expr_str)%+%".\nCheck syntax specification.", call. = FALSE)
  }

  # dprint("evaled_expr: "); dprint(evaled_expr)
  if (!opts$NoChangeFunCalls) {
    # convert one column matrix to a vector:
    f_tovect <- function(X) {
      if (length(dim(X))>1) {
        if (dim(X)[2]==1) X <- as.vector(X)
      }
      X
    }
    if (!is.list(evaled_expr)) {
      evaled_expr <- f_tovect(evaled_expr)
    } else if (is.list(evaled_expr)) {
      evaled_expr <- lapply(evaled_expr, f_tovect)
    }    
  }

  return(list(evaled_expr = evaled_expr, par.nodes = c(Vnames, TD_t_vnames))) # return evaluated expression and parent node names
}

eval.nodeform.asis <- function(expr_call, expr_str, self, data.env) {
  # print("AS IS EVALUTION FOR: "); print(expr_str)
  evaled_expr <- try(eval(expr_call, envir = data.env, enclos = self$user.env)) # eval'ing expr in the envir of data.df
  if(inherits(evaled_expr, "try-error")) {
    stop("error while evaluating node "%+% self$cur.node$name %+%" formula: \n"%+%parse(text = expr_str)%+%".\nCheck syntax specification.", call. = FALSE)
  }
  return(list(evaled_expr = evaled_expr, par.nodes = NULL)) # return evaluated expression and parent node names
}

# ------------------------------------------------------------------------------------------
# **** MOVED THE ENTIRE THING TO R6 CLASS STRUCTURE:
# ------------------------------------------------------------------------------------------
# Function takes a string node formula, current node and current observed data environment
# 1) processes expression into R call, replaces t to its current value
# 2) finds all time-dep var names (Var[]) and non-time dep var names (Var)
# 3) replaces all summary function calls, s.a., func(Var) with apply(Var, 1, func)
# 4) replaces all calls to functions with several vectors, s.a., func(X1,X2,X3) with func(cbind(X1,X2,X3))
# 5) evaluates final expression in a special environment where: 
  # -) variables that have been simulated so far in obs.df are accessible
  # -) the subset vector function '[' is replaces with its specialized version, with syntax TDVar[t_range] for subsetting columns of the observed data by time
  # -) vecapply() function that is a wrapper for apply, converts vector to a 1 col matrix
# 6) standardizes the final expression to be a vector (NEED TO CHANGE FOR CATEGORICAL NODES - sapply over each prob formula in expression?)
eval.nodeform.out <- function(expr.idx, self, data.df) {
  expr_str <- self$exprs_list[[expr.idx]]
  # sVar.name <- self$sVar.expr.names[expr.idx]
  misXreplace <- self$sVar.misXreplace[expr.idx]
  eval.asis <- self$asis.flags[[expr.idx]]

  if (is.character(expr_str) || is.numeric(expr_str)) {
    expr_call <- try(parse(text=expr_str)[[1]])   # parse expression into a call
    if(inherits(expr_call, "try-error")) {
      stop("error while evaluating node "%+% self$cur.node$name %+%" formula:\n "%+%parse(text=expr_str)%+%".\nCheck syntax specification.", call.=FALSE)
    }

  } else if (is.call(expr_str)){
    expr_call <- expr_str
    warning("node "%+%self$cur.node$name%+%": formula is already a parsed call")
  } else {
    stop("node "%+%self$cur.node$name%+%": currently can't process node formulas that are not strings or calls")
  }
  
  # Replace t in the node formula expression with current t value; Replace Kmax its val (returns a call):
  expr_call <- eval(substitute(substitute(e, list(t = eval(self$cur.node$t), Kmax = eval(self$netind_cl$Kmax))), list(e = expr_call)))

  # Removed self$node_fun from data.env as they interfere with R expressions parsing in nodeform_parsers:
  eval.sVar.params <- c(list(self = self),
                        self$df.names(data.df), # special var "ANCHOR_ALLVARNMS_VECTOR_0" with names of already simulated vars
                        list(t = self$cur.node$t), 
                        list(misXreplace = misXreplace), # replacement value for missing network covars
                        list(netind_cl = self$netind_cl),
                        list(nF = self$netind_cl$nF)
                        )
  data.env <- c(eval.sVar.params, data.df)

  if (eval.asis) {
    return(eval.nodeform.asis(expr_call = expr_call, expr_str = expr_str, self = self, data.env = data.env))
  } else {
    return(eval.nodeform.full(expr_call = expr_call, expr_str = expr_str, self = self, data.env = data.env))
  }
}

## ---------------------------------------------------------------------
#' R6 class for parsing and evaluating node R expressions.
#'
#' This \pkg{R6} class will parse and evaluate (in the environment of the input data) the node formulas defined by function
#'  \code{\link[simcausal]{node}}.
#'  The node formula expressions (stored in \code{exprs_list}) are evaluated in the environment of the input data.frame.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{exprs_list}} - Deparsed list of node formula expressions (as strings).
#' \item{\code{user.env}} - Captured user-environment from calls to \code{node} that will be used as enclosing environment during evaluation.
#' \item{\code{cur.node}} - Current evaluation node (set by \code{self$eval.nodeforms()})
#' \item{\code{asis.flags}} - List of flags, \code{TRUE} for "as is" node expression evaluation
#' \item{\code{ReplMisVal0}} - A logical vector that captures args \code{replaceNAw0=TRUE/FALSE} in \code{node} function call.
#'  If \code{TRUE} for a particular node formula in \code{exprs_list} then all missing network \code{VarNode} 
#'  values (when \code{nF[i] < Kmax}) will get replaced with with corresponding value in code{sVar.misXreplace} (default is 0).
#' \item{\code{sVar.misXreplace}} - Replacement values for missing sVar, vector of \code{length(exprs_list)}.
#' \item{\code{netind_cl}} - Pointer to a network instance of class \code{simcausal::NetIndClass}.
#' \item{\code{Kmax}} - Maximum number of friends for any observation.
#' \item{\code{Nsamp}} - Sample size (nrows) of the simulation dataset.
#' \item{\code{node_fun}} - List that contains special subsetting functions \code{'['} and \code{'[['}, where \code{'['} 
#'  is used for subsetting time-varyng nodes and \code{'[['} is used for subsetting network covariate values.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(netind_cl}}{Instantiates new object of class \code{Define_sVar}. 
#'    \code{netind_cl} is the input network stored in an object of class \code{\link[simcausal]{NetIndClass}}.}
#'   \item{\code{set.new.exprs(exprs_list)}}{Sets the internal node formula expressions to the list provided in \code{exprs_list}.}
#'   \item{\code{eval.nodeforms(cur.node, data.df)}}{Evaluate the expressions one by one, returning a list with evaluated expressions. 
#'    \code{cur.node} is the current node object defined with function \code{node} and \code{data.df} is the input data.frame.}
#'   \item{\code{df.names(data.df)}}{List of variables in the input data \code{data.df} gets assigned to a special variable 
#'   (\code{ANCHOR_ALLVARNMS_VECTOR_0}).}
#' }
#' @importFrom assertthat assert_that
#' @export
Define_sVar <- R6Class("Define_sVar",
  class = TRUE,
  portable = TRUE,
  public = list(
    exprs_list = list(),          # expressions converted to strings in a list where attribute "names" is set to the user-supplied vector of expression argument names
    user.env = NULL,              # user environment to be used as enclos arg to eval(sVar)
    cur.node = list(),            # current evaluation node (set by self$eval.nodeforms())
    asis.flags = list(),          # list of flags, TRUE for "as is" node expression evaluation
    ReplMisVal0 = FALSE,          # vector of indicators, for each TRUE sVar.expr[[idx]] will replace all NAs with gvars$misXreplace (0)
    sVar.misXreplace = NULL,      # replacement values for missing sVar, vector of length(exprs_list)
    netind_cl = NULL,
    Kmax = NULL,
    Nsamp = NULL,				          # sample size (nrows) of the simulation dataset

    node_fun = list(
      vecapply = function(X, idx, func) { # custom wrapper for apply that turns a vector X into one column matrix
        if (is.vector(X)) dim(X) <- c(length(X), 1) # returns TRUE only if the object is a vector with no attributes apart from names
        # if (is.atomic(x) || is.list(x)) dim(X) <- c(length(X), 1) # alternative way to test for vectors
          x <- parse(text = deparse(func))[[1]]
          nargs <- length(x[[2]])
          if (nargs>1) {
            funline <- deparse(func)[1]
            stop(funline%+%
            ". Node formulas cannot call non-vectorized functions with more than one named argument. If this is a vectorized function, pass its name to set.DAG(, vecfun=).")
          }
        apply(X, idx, func)
      },

      cbind_mod = function(...) { # cbind wrapper for c(,) calls in node formulas, turns one row matrix into repeat Nsamp row matrix
        env <- parent.frame()
        cbind_res <- do.call("cbind", eval(substitute(alist(...)), envir = env) , envir = env)
        if (nrow(cbind_res)==1) {
          # Nsamp <- get("Nsamp", envir = env)
          Nsamp <- env$self$Nsamp
          dprint("env$self$Nsamp:"); dprint(env$self$Nsamp)
          assert_that(!is.null(Nsamp))

          cbind_res <- matrix(cbind_res, nrow = Nsamp, ncol = ncol(cbind_res), byrow = TRUE)
        }
        dprint("cbind_res"); dprint(cbind_res)
        cbind_res
      },

      # custom function for vector look up '['
      # function takes the name of the TD var and index vector => creates a vector of time-varying column names in df
      # returns matrix TD_var[indx]
      # ***NOTE: current '[' cannot evalute subsetting that is based on values of other covariates such as A1C[ifelse(BMI<5, 1, 2)]
      `[` = function(var, indx, ...) {
        env <- parent.frame()
        t <- env$t # t <- get("t", envir = env)
        var <- substitute(var)
        var.chr <- as.character(var)

        if (missing(indx)) stop("missing tindex when using Var[tindex] inside the node formula")
        if (identical(class(indx),"logical")) indx <- which(indx)
        if (is.null(t)) stop("references, s.a. Var[t] are not allowed when t is undefined")
        if (max(indx)>t) stop(paste0(var, "[", max(indx),"] cannot be referenced in node formulas at t = ", t))  # check indx<= t

        # ******* NOTE *******
        # Don't like the current implementation that defines TDvars as characters and then returns a matrix by cbinding 
        # the existing columins in existing data.frame. This is possibly wasteful. Could we instead subset the existing data.frame?
        TDvars <- var.chr%+%"_"%+%indx
        # Checking the variables paste0(var, "_", indx) exist in simulated data.frame environment:
        dprint("ANCHOR_ALLVARNMS_VECTOR_0:"); dprint(env[["ANCHOR_ALLVARNMS_VECTOR_0"]])

        # TO DO: **** 
        # EXTEND TO CHECKING FOR TDvar IN ENCLOSING ENVIRONMENT (user.env) AS WELL IF TDvar_t doesn't exist in the data
        # IF TDvar exists check that its a vector of appropriate length, index it accordinly (using which(t%in%tvec))
        # will need to first eval such vector the variable as in: 
        # var.val <- eval(var, envir = env)
        existsTDVar <- function(TDvar_t) TDvar_t %in% env[["ANCHOR_ALLVARNMS_VECTOR_0"]]
        check_exist <- sapply(TDvars, existsTDVar)
        if (!all(check_exist)) stop("undefined time-dependent variable(s): "%+%TDvars[which(!check_exist)])
        # THIS STEP COULD BE MORE MEMORY EFFICIENT IF WAS SUBSETTING INSTEAD (BY COLS) ON EXISTING data MATRIX:
        TDvars_eval <- eval(parse(text=paste0("cbind(",paste0(TDvars, collapse=","),")")), envir = env)
        return(TDvars_eval)
      },

      # Builds netVar matrix by using matrix env$NetIndobj$NetInd_k, cbind on result
      # For W[[0]] to work without if else below need to do this:
      # NetInd_k <- cbind(c(1:n), NetInd_k) and then netidx <- netidx + 1
      `[[` = function(var, netidx, ...) {
        env <- parent.frame()
        t <- env$t # t <- get("t", envir = env)
        if (!is.null(t)) stop("simultaneous time varying node references Var[t] and network references Var[[netidx]] are currently not supported")
        if (missing(netidx)) stop("network index (netidx) must be specified when using Var[[netidx]]")
        netind_cl <- env$netind_cl
        if (is.null(netind_cl)) stop("Network must be defined when using Var[[netidx]] syntax")
        Kmax <- netind_cl$Kmax

        var <- substitute(var)
        var.chr <- as.character(var)
        if (! (var.chr %in% env[["ANCHOR_ALLVARNMS_VECTOR_0"]])) stop("variable " %+% var.chr %+% " doesn't exist")
        var.val <- eval(var, envir = env)
        n <- length(var.val)
        if (identical(class(netidx),"logical")) netidx <- which(netidx)
        netVars_eval <- matrix(0L, nrow = n, ncol = length(netidx))
        colnames(netVars_eval) <- netvar(var.chr, netidx)
        for (neti in seq_along(netidx)) {
          if (netidx[neti] %in% 0L) {
            netVars_eval[, neti] <- var.val
          } else {
            netVars_eval[, neti] <- var.val[netind_cl$NetInd_k[, netidx[neti]]]
            # opting for replace on entire netVars_eval, will need to do benchmarks later to compare:
            # netVars_eval[is.na(netVars_eval[, neti]), neti] <- env$misXreplace
          }
        }
        # Don't need to do this if env$misXreplace==gvars$misval (i.e., when want to leave NAs as is)
        netVars_eval[is.na(netVars_eval)] <- env$misXreplace
        return(netVars_eval)
      }
    ),

    initialize = function(netind_cl) {
      self$netind_cl <- netind_cl
      self$Kmax <- self$netind_cl$Kmax
      invisible(self)
    },

    set.new.exprs = function(exprs_list) {
      self$exprs_list <- exprs_list
      self$asis.flags <- attributes(exprs_list)[["asis.flags"]]
      # check for special argument replaceNAw0, if exists, remove it from the list of expressions:
      if (any(names(self$exprs_list) %in% "replaceNAw0")) {
        ReplMisVal0.idx <- which(names(self$exprs_list) %in% "replaceNAw0")
        self$ReplMisVal0 <- as.logical(self$exprs_list[[ReplMisVal0.idx]])
        self$exprs_list <- self$exprs_list[-ReplMisVal0.idx]
        print("Detected replaceNAw0 flag with value: " %+% self$ReplMisVal0);
      }
      # if doesn't already exist, init setting for the names attribute of self$exprs_list:
      if (is.null(names(self$exprs_list))) names(self$exprs_list) <- rep_len("", length(self$exprs_list))
      # OS 09/22/15: Disabling the error for un-named expressions (name(self$exprs_list[[i]])=="");
      # (1) Want to allow un-named expressions in tmlenet::def.sW & tmlenet::def.sA;
      # (2) This condition is already pre-checked inside simcausal::node;
      # if (length(self$exprs_list) != 0 && (is.null(self$sVar.expr.names) || any(self$sVar.expr.names==""))) {
      #   stop("must provide a name for each node expression")
      # }

      if (is.null(self$asis.flags)) { 
        self$asis.flags <- as.list(rep.int(FALSE, length(self$exprs_list)))
        names(self$asis.flags) <- names(self$exprs_list)
      }
      self$ReplMisVal0 <- rep_len(self$ReplMisVal0, length(self$exprs_list))
      self$sVar.misXreplace <- ifelse(self$ReplMisVal0, gvars$misXreplace, gvars$misval)
      # self$sVar.noname <- rep_len(self$sVar.noname, length(self$exprs_list))
      print("Final node expression(s) list: "); print(self$exprs_list)
      invisible(self)
    },

    eval.nodeforms = function(cur.node, data.df) {
      assert_that(is.environment(self$user.env))
      self$cur.node <- cur.node
      self$ReplMisVal0 <- FALSE
      # Parse the node formulas (parameters), set self$exprs_list, set new self$sVar.misXreplace and self$sVar.noname if found:
      self$set.new.exprs(exprs_list = cur.node$dist_params)
      self$Nsamp <- nrow(data.df)
      sVar.res_l <- lapply(seq_along(self$exprs_list), eval.nodeform.out, self = self, data.df = data.df)
      names(sVar.res_l) <- names(self$exprs_list)
      return(sVar.res_l)
    },

    df.names = function(data.df) { # list of variable names from data.df with special var name (ANCHOR_ALLVARNMS_VECTOR_0)
      return(list(ANCHOR_ALLVARNMS_VECTOR_0 = colnames(data.df)))
    },

    # This user.env is used for eval'ing each sVar exprs (enclos = user.env)
    set.user.env = function(user.env) {
      assert_that(!is.null(user.env))
      assert_that(is.environment(user.env))
      self$user.env <- user.env
    }
  ),

  active = list(
    placeholder = function() {}
  ),

  private = list(
    privplaceholder = function() {}
  )
)
# nocov end