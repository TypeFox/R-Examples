# This file contains functions for testing:
#       do.test
#       allTrue
#       expectStop
#       expectWarnings
# See also:
#       help(do.test) contains all.equal.excluding

do.test <- 
function(file, verbose = FALSE, strict = FALSE, local = FALSE, check)
{
  # This is a port of the S-PLUS version of do.test to R

  # Expressions will be parsed and evaluated from `file'.
  # Each expression should evaluate to a logical TRUE.
  # Otherwise, `do.test()' will print the expression and its
  # value.  If `verbose' is `TRUE', the printing is done anyway.
  # If `strict' is `TRUE', any validity failures cause an error;
  # that is, you get to debug after the first failed assertion.
  # If `check' is supplied, `do.test' evaluates this expression
  # (it should be given via `Quote()' between each parse and
  # evaluation.  (This is for when you need to check some global
  # information.)  The `local' argument controls where the
  # evaluation takes place: by default, in the caller of do.test,
  # typically the global environment (objects created remain there
  # after `do.test' is finished).
  # `local=T' causes `do.test' to create and work in a new environment.
  if(is(file, "connection") || (is.R() && inherits(file, "connection")))
    # R's is() currently ignores inheritance
    cat("---- Test connection:",
	if(is.R()) summary(file)$description else file@description, 
	"Class:", class(file), "----\n") else {
    cat("------- Test file:", file, "---------\n")
  }
  E <- parse(file = file, n = -1)
  # For R, printing expressions that fail would miss comments;
  # comments are stored with source text in an attribute "srcref".
  if(is.R())
    E2 <- attr(E, "srcref")

  Env <- if(local) new.env() else .GlobalEnv
  for(I in seq(along = E)) {
    if(!missing(check))
      eval(check, envir=Env) #eval(check, local)
    if(verbose)
      print(if(is.R()) E2[[I]] else E[[I]])
    val <- try(eval(E[[I]], envir=Env), silent=TRUE)
    failed <- (!identical(val, TRUE))
    if(verbose || failed) {
      if(!verbose)
	print(if(is.R()) E2[[I]] else E[[I]])
      print(val)
      if(strict && failed)
        stop("Test Failed")
    }
  }
  invisible(NULL)
}
# Use silent=TRUE to match S+ default
# If an expression fails, print E2 (from srcref attribute) to keep comments.

allTrue <- 
function(...)
{
  # arguments are typically results from all.equal()
  # return TRUE if all are TRUE; else show results which are not true
  #
  # This function is useful in loop tests, to do multiple comparisons
  # { object1 <- (some code)
  #   object2 <- (other code)
  #   allTrue(all.equal(object1$a, object2$a), 
  #      all.equal(object1$b, object2$b), 
  #      all.equal(object1$c, object2$c))
  # }
  # Then if any of the tests fail, the whole loop test is printed, 
  # including the calls that created the data.
  L <- list(...)
  il <- sapply(L, function(x) is.logical(x) && !is.na(x) && x)
  if(all(il))
    return(TRUE)
  # At least one test fails; indicate which tests fail.
  if(is.null(names(L)))
    names(L) <- paste("Comparison", seq(length = length(L))
      ) else for(i in 1:length(L)) {
	if(is.null(names(L)[i]))
	  names(L)[i] <- paste("Comparison", i)
      }
  L[!il]
}

expectStop <- 
function(expr, expected = NULL)
{
  # expr is an expression (unquote) to be evaluated, typically a compound
  #   expression like those typically found in loop test files called by
  #   do.test.  
  # expected is either NULL, or a text string containing part of the
  #   expected message.  This may be a regular expression in the
  #   format used by regmatch.
  #
  # This function is useful for checking error checking; that a function
  # stops when it should, and gives the right message.
  # For example, this may be in a file called by do.test:
  #     {
  #       expectStop(var(1:5, 1:4), 
  #                  "x and y must have the same number of")
  #     }
  #   
  # The function returns TRUE if 
  #   1) a stop() occurs
  #   2) the error message is expected
  # Otherwise it returns appropriate messages.
  expr <- substitute(msg <- try(EXPR, silent=TRUE),
		     list(EXPR = substitute(expr)))
  msg <- eval(expr, sys.parent())
  if(!is(msg, "try-error"))
    return(paste("Failed to stop, instead returned an object of class",
		 paste(class(msg), collapse=", ")))
  if(is.null(expected))
    return(TRUE)
  if(length(grep(expected, msg[1])))
    return(TRUE)
  c("Error message differs from expected message", "Expected:", expected, 
    "Actual:", msg[1])
}

expectWarnings <- function(expr, expected){
  # Evaluate expr, and compare the actual warnings to expected warnings
  #
  # expr is an expression (unquoted) to be evaluated, typically a compound
  #   expression like those typically found in loop test files called by
  #   do.test.
  # expected is a vector of expected warnings (text strings, may
  #   be abbreviated)
  #
  # This function is useful for detecting expected warning messages in
  # loop tests.  For instance, the warnings generated by the loop test
  # expression
  #
  # {
  #   object1 <- (code generating warning messages)
  #   object2 <- (code generating possibly other warning messages)
  #   all.equal(object1, object2)
  # }
  #
  # can be verified and suppressed by using
  #
  # {
  #   expectWarnings(
  #   {
  #     object1 <- (code generating warning messages)
  #     object2 <- (code generating possibly other warning messages)
  #     all.equal(object1, object2)
  #   },
  #   c("expected warning 1",
  #     "expected warning 2",
  #     ...)
  # }
  #
  # The function returns TRUE if
  #   1) expr evaluates to TRUE; and
  #   2) each warning message produced by evaluating expr contains as a
  #      substring an element of expected, and each element of expected
  #      is a substring of at least one of the produced warning messages.
  #
  # Otherwise it returns a list with some or all of the following components:
  #
  #  "Test result" -- the value (if not TRUE) returned by expr.
  #
  #  "Unexpected warnings" -- generated warning messages not found in
  #      expected.
  #
  #  "Warnings expected but not found" -- messages in expected not generated
  #      by evaluating expr.
  #
  # Normal printing of warning messages is suppressed.
  #
  # evaluate expr,
  # check the actual warnings list against the expected
  # return differences
  ENV <- environment()
  actualWarnings <- NULL
  warningsHandler <- function(w) {
    assign("actualWarnings",
	   c(get("actualWarnings", envir=ENV), w$message),
	   envir = ENV)
    invokeRestart("muffleWarning")
  }
  value <- withCallingHandlers(eval(substitute(expr), envir = parent.frame()),
                    warning = warningsHandler)

  # Return a combination of value and a summary of the discrepancies
  # between actual and expected warnings
  warnings2 <- unique(actualWarnings)
  if(lw2 <- length(warnings2)) {
    matches <- charmatch(expected, warnings2)
    w <- is.na(matches)
    warnNotFound <- expected[w]
    # handle ambiguous matches individually
    w0 <- matches == 0
    matches <- matches[!(w | w0)]
    for(i in which(w0)) {
      for(j in 1:lw2) {
        if(pmatch(expected[i], warnings2[j], nomatch = FALSE))
          matches <- c(matches, j)
      }
    }
    w2 <- setdiff(seq(along = warnings2), unique(matches))
    warnUnexpected <- warnings2[w2]
  } else {
    warnNotFound <- expected
    warnUnexpected <- NULL
  }
  allTrue("Test result" = value,
	  "Unexpected warnings" =
	  if(length(warnUnexpected)) warnUnexpected else TRUE,
	  "Warnings expected but not found" =
	  if(length(warnNotFound)) warnNotFound else TRUE)
}
