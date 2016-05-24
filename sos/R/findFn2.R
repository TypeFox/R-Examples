`?` <- function (e1, e2)
{
  call <- match.call()

  original <- function() {
  # call the original ? function
    call[[1]] <- quote(utils::`?`)
    return(eval(call, parent.frame(2)))
  }

  # We don't handle requests with type
  if (!missing(e2)) {
    return(original())
  }

# We only handle function calls with double ??
#    (not counting the original one)
  topicExpr1 <- substitute(e1)
  if (!is.call(topicExpr1)
      || length(topicExpr1) != 2
      || topicExpr1[[1]] != "?"
      || !is.call(topicExpr1[[2]])
      || length(topicExpr1[[2]]) != 2
      || topicExpr1[[2]][[1]] != "?")
    return(original())

 # Get the expression
  topicExpr <- topicExpr1[[2]][[2]]


  # Construct our call to RSiteSearch.function
  if (is.call(topicExpr)) {
# It must not be a call to ?,
#      that would mean there are 4 or more
    if (topicExpr[[1]] == "?")return(original())
    lastArg <- length(topicExpr)
    topicExpr[[lastArg+1]] <- as.character(topicExpr[[1]])
    names(topicExpr)[[lastArg+1]] <- "string"
#    topicExpr[[1]] <- quote(sos::findFn)
    topicExpr[[1]] <- quote(findFn)
    f. <- eval(topicExpr, parent.frame(1))
  } else {
#    	RSiteSearch.function(as.character(topicExpr))
#    ff <- findFn(as.character(topicExpr))
    f. <- do.call('findFn', list(as.character(topicExpr)))
  }
  f.
}
