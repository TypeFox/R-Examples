# Adds post-assertions to the base function
'%must%' <- function(fn.ref, condition)
{
  strict <- TRUE
  child <- deparse(substitute(fn.ref))

  expr <- deparse(substitute(condition))
  if (length(grep('function', expr)) < 1) 
  {
    # TODO: Version 1 is deprecated
    if (paradigm.options('version') == 1)
      return(.ensure(child, expr, strict, label='ensure.xps'))

    # This is version 2
    return(.must(child, expr))
  }

  return(.ensure(child, condition, strict, label='ensure.fns'))
}

.must <- function(parent, condition)
{
  # We use 2 because this is called from within the 'guard' function so the
  # stack is two down
  where <- topenv(parent.frame(2))
  .setup.parent(parent, where)

  fn <- get(parent, where)
  variant.count <- attr(fn,'variant.count')
  name <- paste(parent, variant.count, collapse=".")
  gs <- attr(fn, name)
  gs$ensures <- c(gs$ensures, condition)

  attr(fn, name) <- gs
  assign(parent, fn, where)
  invisible()
}

# Shortcut form for expressions instead of more verbose functions
# This is the standard way of writing ensures, although the long form is
# slightly more efficient from an execution perspective
# label := { ensure.xps, ensure.fns }
.ensure <- function(child, condition, strict, label)
{
  # Need to add a lazy binding so that the function is bound only when
  # the actual child.fn is defined. Otherwise, it's impossible to determine
  # the arguments
  
  parent <- sub('\\.[^.]+$','', child)
  where <- paradigm.options(parent)
  # We use 2 because this is called from within the 'ensure' function so the
  # stack is two down
  if (is.null(where)) where <- topenv(parent.frame(2))

  if (! exists(parent, where))
  {
    msg <- "Function %s has no visible parent function '%s'"
    stop(sprintf(msg, child, parent))
  }
  fn <- get(parent, where)
  gs <- attr(fn, label)
  if (is.null(gs)) gs <- list()
  #gs[[child]] <- c(gs[[child]], expr)
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

# Unlike guards, ensures always execute following a matching guard, so no
# sophisticated logic is needed here.
.applyEnsure <- function(ensures, child, validate.fn, result, ...)
{
  if (is.null(ensures)) return(0)

  args <- list(...)
  idx <- 0
  valid <- TRUE
  f <- child
  for (g in ensures[[child]])
  {
    idx <- idx + 1
    valid <- valid && validate.fn(g,f, result, ...)

    if (! valid) break
  }
  ifelse(valid, 0, idx)
}

# This only works for one level
ensures <- function(fn, inherits=TRUE, child=NULL)
{
  if (! is.function(fn))
    stop("Ensure introspection can only be applied to functions")

  gfs <- attr(fn, 'ensure.fns', exact=TRUE)
  gxs <- attr(fn, 'ensure.xps', exact=TRUE)
  gs <- list(functions=gfs, expressions=gxs)

  if (! is.null(gfs) || ! is.null(gxs)) return(gs)
  if (! inherits) return(gs)

  parent <- sub('\\.[^.]+$','', deparse(substitute(fn)) )

  # At root of hierarchy
  if (length(grep('^get\\(', parent)) > 0) return(NULL)

  ensures(get(parent, inherits=TRUE), child=child)
}


.validateEnsureFunction <- function(g,f, result, ...)
{
  if (! is.function(g)) stop("Invalid ensure function specified")
  g(result, ...)
}

.ENSURE_EXPRESSION <- "%s <- function(result,%s) { %s }"
.validateEnsureExpression <- function(g,f, result, ...)
{
  f.exec <- get(f)
  fn.handle <- paste('.ensure',f,sep='_')
  my.args <- paste(names(formals(f.exec)), collapse=',')
  xps <- parse(text=sprintf(.ENSURE_EXPRESSION, fn.handle,my.args, g))
  eval(xps)
  eval(parse(text=sprintf("%s(result, ...)",fn.handle)))
}

