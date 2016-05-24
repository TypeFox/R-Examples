#library(futile.logger)
#configLogger(threshold=DEBUG)
#lg <- getLogger()
# Probably don't need this any more
paradigm.options <- OptionsManager('paradigm.options', list(version=2))


# Adds guards to the base function for functional dispatching
'%when%' <- function(fn.ref, condition)
{
  strict <- TRUE
  child <- deparse(substitute(fn.ref))

  expr <- deparse(substitute(condition))
  if (length(grep('^(c\\()?function', expr, perl=TRUE)) < 1) 
  {
    # TODO: Version 1 is deprecated
    if (paradigm.options('version') == 1)
      return(.guard(child, expr, strict, label='guard.xps'))

    # This is version 2
    return(.when(child, expr))
  }

  return(.guard(child, condition, strict, label='guard.fns'))
}

# Adds guards to the same function variant. This is only for version 2
'%also%' <- function(fn.ref, condition)
{
  child <- deparse(substitute(fn.ref))
  expr <- deparse(substitute(condition))
  return(.also(child, expr))
}

# Adds the function definition to the function
'%as%' <- function(fn.ref, fn.def)
{
  child <- deparse(substitute(fn.ref))
  return(.as(child, fn.def))
}

'%default%' <- function(fn.ref, fn.def)
{
  child <- deparse(substitute(fn.ref))
  return(.default(child, fn.def))
}


.setup.parent <- function(parent, where)
{
  # Overwrite a final definition (as opposed to appending)
  if (exists(parent, where))
  {
    parent.def <- get(parent, where)
    is.final <- attr(parent.def, 'final')
    if (!is.null(is.final) && is.final == TRUE)
    {
      attributes(fn) <- NULL
      attr(parent.def,'variant.count') <- 0
      assign(parent, parent.def, where)
    }
  }
  else if (! exists(parent, where))
  {
    pattern <- 'function(...) UseFunction(\'%s\',...)'
    parent.def <- eval(parse(text=sprintf(pattern,parent)))
    attr(parent.def,'variant.count') <- 0
    #cat("Adding parent function",parent.def,"to",where,"\n")
    assign(parent, parent.def, where)
    #msg <- "Function %s has no visible parent function '%s'"
    #stop(sprintf(msg, child, parent))
  }
}

# New structure is
# attr(,"variant.count")
# attr(,"parent.1")$guards
# attr(,"parent.1")$ensures
# attr(,"parent.1")$definition
# attr(,"parent.2")$guards
# attr(,"parent.2")$definition
.when <- function(parent, condition)
{
  # We use 2 because this is called from within the 'guard' function so the
  # stack is two down
  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)

  fn <- get(parent, where)
  variant.count <- attr(fn,'variant.count') + 1
  name <- paste(parent, variant.count, sep=".")
  #cat("[.when] Adding",name,"to function\n")
  gs <- list(guards=c(condition))

  attr(fn, name) <- gs
  attr(fn, 'variant.count') <- variant.count
  assign(parent, fn, where)
  invisible()
}

.also <- function(parent, condition)
{
  # We use 2 because this is called from within the 'guard' function so the
  # stack is two down
  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)

  fn <- get(parent, where)
  variant.count <- attr(fn,'variant.count')
  name <- paste(parent, variant.count, sep=".")
  gs <- attr(fn, name)
  gs$guards <- c(gs$guards,condition)

  attr(fn, name) <- gs
  assign(parent, fn, where)
  invisible()
}

.as <- function(parent, fn.def)
{
  # We use 2 because this is called from within the 'guard' function so the
  # stack is two down
  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)

  fn <- get(parent, where)
  variant.count <- attr(fn,'variant.count')
  name <- paste(parent, variant.count, sep=".")
  gs <- attr(fn, name)
  count <- length(gs)
  gs$definition <- fn.def

  attr(fn, name) <- gs
  assign(parent, fn, where)
  invisible()
}

.default <- function(parent, fn.def)
{
  # We use 2 because this is called from within the 'guard' function so the
  # stack is two down
  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)

  fn <- get(parent, where)
  attr(fn, 'default') <- fn.def
  assign(parent, fn, where)
  invisible()
}


# DEPRECATED
.guard <- function(child, condition, strict, label)
{
  parent <- sub('\\.[^.]+$','', child)
  where <- paradigm.options(parent)
  # For debugging
  #  cat("[From .guard]\n")
  #  cat("Parent: "); print(parent.frame(2))
  #  cat("TopEnv: "); print(topenv(parent.frame(2)))

  # We use 2 because this is called from within the 'guard' function so the
  # stack is two down
  if (is.null(where)) where <- topenv(parent.frame(2))
  # This is only the frame and not the environment
  #if (is.null(where)) where <- parent.frame(2)
  # This doesn't work
  #if (is.null(where)) where <- topenv(parent.frame())

  if (! exists(parent, where))
  {
    pattern <- 'function(...) UseFunction(\'%s\',...)'
    parent.def <- eval(parse(text=sprintf(pattern,parent)))
    #cat("Adding parent function",parent.def,"to",where,"\n")
    assign(parent, parent.def, where)
    #msg <- "Function %s has no visible parent function '%s'"
    #stop(sprintf(msg, child, parent))
  }
  fn <- get(parent, where)
  gs <- attr(fn, label)
  if (is.null(gs)) gs <- list()
  #gs[[child]] <- c(gs[[child]], condition)
  gs[[child]] <- c(condition)
  attr(fn, label) <- gs

  if (strict)
  {
    ss <- attr(fn, 'strict')
    if (is.null(ss)) ss <- list()
    ss[[child]] <- strict
    attr(fn, 'strict') <- ss
  }

  assign(parent, fn, where)

  invisible()
}

# Get the guards for a function. The function can be either the parent or any
# of its attached children
# This only works for one level since the deparse/substitute's lazy evaluation
# messes up the recursion
guards <- function(fn.ref, inherits=TRUE)
{
  if (! is.function(fn.ref))
    stop("Guard introspection can only be applied to functions")

  if (paradigm.options('version') == 1)
  {
    gfs <- attr(fn.ref, 'guard.fns', exact=TRUE)
    gxs <- attr(fn.ref, 'guard.xps', exact=TRUE)
    gs <- list(functions=gfs, expressions=gxs)
    if (! is.null(gfs) || ! is.null(gxs)) return(gs)
    if (! inherits) return(gs)

    parent <- sub('\\.[^.]+$','', deparse(substitute(fn.ref)))
    guards(get(parent, inherits=TRUE), inherits=TRUE)
  }
  else
  {
    attributes(fn.ref)
  }
}

# Returns a function variant. Useful for debugging
# Do not enclose the name in quotes or this will fail
# debug(variant(some_function.1))
variant <- function(name.fn)
{
  name <- deparse(substitute(name.fn))
  parts <- strsplit(name, ".", fixed=TRUE)[[1]]
  parent <- parts[1:(length(parts)-1)]
  #index <- parts[length(parts)]

  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)
  fn <- get(parent)
  attr(fn, name)$definition
}

rm.variant <- function(name.fn)
{
  name <- deparse(substitute(name.fn))
  parts <- strsplit(name, ".", fixed=TRUE)[[1]]
  parent <- parts[1:(length(parts)-1)]

  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)
  fn <- get(parent)
  attr(fn, name) <- NULL
  assign(parent, fn, where)
}

# Operates on a child function or function name
# DEPRECATED
.isStrict <- function(child.fn)
{
  if (is.function(child.fn)) child <- deparse(substitute(child.fn))
  else child <- child.fn

  parent <- sub('\\.[^.]+$','', child)
  fn <- get(parent, inherits=TRUE)

  ss <- attr(fn, 'strict', exact=TRUE)
  #cat("Got strict on function:",ss[[child]],"\n")
  if (! is.null(ss) && ! is.null(ss[[child]])) return(ss[[child]])

  FALSE
}

'%isa%' <- function(argument, type)
{
  type <- gsub('[\'"]','',deparse(substitute(type)))
  type %in% class(argument)
}

# Note this will produce a vector of results
'%hasa%' <- function(argument, property)
{
  property <- gsub('[\'"]','',deparse(substitute(property)))
  property <- gsub(' ','', property, fixed=TRUE)
  property <- sub('c(','', property, fixed=TRUE)
  property <- sub(')','', property, fixed=TRUE)
  props <- strsplit(property, ',', fixed=TRUE)[[1]]
  props %in% names(argument)
}

'%hasall%' <- function(argument, property)
{
  property <- gsub('[\'"]','',deparse(substitute(property)))
  property <- gsub(' ','', property, fixed=TRUE)
  property <- sub('c(','', property, fixed=TRUE)
  property <- sub(')','', property, fixed=TRUE)
  props <- strsplit(property, ',', fixed=TRUE)[[1]]
  all(props %in% names(argument))
}


.SIMPLE_TYPES <- c('numeric','character','POSIXt','POSIXct')
.is.simple <- function(x) any(class(x) %in% .SIMPLE_TYPES)
.as.simple <- function(x)
{
  if (! .is.simple(x)) return(class(x)[1])
  if (length(x) == 1) return(x)
  if (length(x) < 5) sprintf("c(%s)", paste(x, collapse=','))
  else sprintf("c(%s, ...)", paste(x[1:4], collapse=','))
}

# Dispatcher for a more functional paradigm. This executes a function based on
# which guards are matched. The order of evaluation is based on the order the
# guards are declared.
.ERR_USE_FUNCTION <- "No valid function for '%s/%s' : %s"
.ERR_ENSURE_FAILED <- "Assertion '%s' failed for args = %s and result = %s"
UseFunction <- function(fn.name, ...)
{
  result <- NULL
  #cat("[UseFunction] Getting function",fn.name,"\n")
  fn <- get(fn.name)
  #cat("[UseFunction] Getting guards for",fn.name,"\n")
  gs <- guards(fn, inherits=FALSE)
  if (is.null(gs)) stop("Incorrect guard output. Please report to maintainer.")
  # TODO: Version 1 is deprecated
  if (paradigm.options('version') == 1)
  {
    if (is.null(gs$functions) && is.null(gs$expressions))
      stop("Function must have guards for functional dispatching")

    matched.fn <- .applyGuard(gs$functions, .validateGuardFunction, ...)
    if (is.null(matched.fn))
      matched.fn <- .applyGuard(gs$expressions, .validateGuardExpression, ...)
  }
  else
  {
    #cat("[UseFunction] Looking for matching variant for",fn.name,"\n")
    matched.fn <- .applyVariant(gs, .validateGuardExpression, ...)
  }
  if (is.null(matched.fn))
  {
    if (! is.null(attr(fn, 'default')))
    {
      matched.fn <- attr(fn, 'default')
    }
    else
    {
      args <- sapply(list(...), .as.simple)
      arg.names <- paste(args, collapse=', ')
      arg.length <- length(args)
      stop(sprintf(.ERR_USE_FUNCTION, fn.name, arg.length, arg.names))
    }
  }

  #cat("[UseFunction] Calling matched function\n")
  #print(matched.fn)
  #cat("\n")
  result <- do.call(matched.fn, list(...))

  es <- ensures(fn, inherits=FALSE, child=matched.fn)
  if (is.null(es)) return(result)
  if (is.null(es$functions) && is.null(es$expressions)) return(result)

  idx <- .applyEnsure(es$functions, matched.fn,
    .validateEnsureFunction, result=result, ...)
  if (idx > 0)
  {
    err.msg <- sprintf(.ERR_ENSURE_FAILED, es$functions[idx], 
      paste(sapply(list(...), .as.simple), collapse=', '), .as.simple(result))
    stop(err.msg)
  }

  idx <- .applyEnsure(es$expressions, matched.fn,
    .validateEnsureExpression, result=result, ...)
  if (idx > 0)
  {
    err.msg <- sprintf(.ERR_ENSURE_FAILED, es$expressions[idx], 
      paste(sapply(list(...), .as.simple), collapse=', '), .as.simple(result))
    stop(err.msg)
  }
      

  result
}

.validateGuardFunction <- function(g,f, ...)
{
  if (is.logical(g)) g
  else if (is.function(g)) g(...)
  else NULL
}

.GUARD_EXPRESSION <- "%s <- function(%s) { %s }"
.validateGuardExpression <- function(g,f.exec, ...)
{
  f <- runif(1)
  fn.handle <- paste('.guard',f,sep='_')
  my.args <- paste(names(formals(f.exec)), collapse=',')
  xps <- parse(text=sprintf(.GUARD_EXPRESSION, fn.handle,my.args, g))
  eval(xps)
  #eval(parse(text=sprintf("%s(...)",fn.handle)))
  eval(parse(text=sprintf("%s(...)",fn.handle)))
}

    # Experimental to try to include default arguments
    #if (length(args) > length(formals(f.exec)) ) next
    # This is the number of arguments without default values
    #required <- sapply(formals(f.exec), function(x) x == '')
    #min.args <- sum(ifelse(required,1,0))
    #if (length(args) < min.args) next

.applyVariant <- function(defs, validator, ...)
{
  args <- list(...)
  for (f in names(defs))
  {
    #cat("[.applyVariant] Checking variant",f,"\n")
    if (class(defs[[f]]) != 'list') next
    #cat("[.applyVariant] Elements:\n")
    #print(defs[[f]])

    #cat("[.applyVariant] Getting function reference",f,"\n")
    f.exec <- defs[[f]][['definition']]
    # Basic validation
    if (is.null(f.exec)) next
    if (length(formals(f.exec)) != length(args)) next

    # Matched named arguments
    #cat("[.applyVariant] Getting arguments for",f,"\n")
    non.empty <- names(args)[nchar(names(args)) > 0]
    #cat("[.applyVariant] Matching arguments for",f,"\n")
    if (length(non.empty) > 0 &&
        length(setdiff(non.empty, names(formals(f.exec)))) > 0 ) next

    valid <- TRUE
    #cat("[.applyVariant] Checking guards for",f,"\n")
    for (g in defs[[f]][['guards']])
    {
      #cat("[.applyVariant] Checking guard",g,"\n")
      this.valid <- validator(g,f.exec, ...)
      #cat("[.applyVariant] Got result",this.valid,"\n")
      if (is.null(this.valid))
      {
        msg <- "Skipping invalid guard '%s' for function '%s'"
        cat(sprintf(msg, g, f))
        next
      }

      valid <- valid && this.valid
      # NOTE: If later on we want guard sequences instead, test for valid
      # instead so the logic becomes a disjunction
      if (! valid) break
    }
    if (valid) return(f.exec)
  }
  NULL
}

.applyGuard <- function(guards, validator, ...)
{
  args <- list(...)
  # First iterate through formal functions. Formal functions have precedence
  # over expressions.
  for (f in names(guards))
  {
    # First check that this function is valid and matches the argument count
    f.exec <- get(f)
    if (is.null(f.exec)) next
    if (length(formals(f.exec)) != length(args)) next

    # If strict, match exactly the function arguments with the arguments
    # passed in. This is the default behavior.
    non.empty <- names(args)[nchar(names(args)) > 0]
    if (length(non.empty) > 0 && .isStrict(f) && 
        length(setdiff(non.empty, names(formals(f)))) > 0 )
      next

    valid <- TRUE
    for (g in guards[[f]])
    {
      this.valid <- validator(g,f.exec, ...)
      if (is.null(this.valid))
      {
        msg <- "Skipping invalid guard '%s' for function '%s'"
        cat(sprintf(msg, g, f))
        next
      }

      valid <- valid && this.valid
      # NOTE: If later on we want guard sequences instead, test for valid
      # instead so the logic becomes a disjunction
      if (! valid) break
    }
    if (valid) return(f)
  }
  NULL
}

# General framework for accessing an object. For create and access, strings can
# be used to identify the class type. These don't use UseMethod for syntactic
# sugar inconsistencies (UseMethod passes original value as opposed to 
# transformed values (i.e. from symbol to character).
# The idea for AbuseMethod is to define high-level framework APIs that can be 
# extended by user code. This produces some consistency in the way common
# operations are called. Functions dispatched from AbuseMethod can contain 
# UseFunction declarations for further dispatching.
AbuseMethod <- function(fn.name, type, ..., EXPLICIT=FALSE, ALWAYS=TRUE)
{
  if (EXPLICIT)
  {
    target <- paste(fn.name,type,sep='.')
    if (! exists(target))
    {
      if (ALWAYS) target <- paste(fn.name,'default',sep='.')
      else stop("No valid method found")
    }

    do.call(target, list(...) )
  }
  else
  {
    types <- class(type)
    for (t in types)
    {
      target <- paste(fn.name,t,sep='.')
      if (exists(target)) break
    }
    if (! exists(target))
    {
      if (ALWAYS) target <- paste(fn.name,'default',sep='.')
      else stop("No valid method found")
    }
    do.call(target, c(type, list(...)) )
  }
}



