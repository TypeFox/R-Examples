# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

.resampleMakeFunction <- function(data, statistic,
                                  substituteData, substituteStatistic,
                                  args.stat = NULL,
                                  resampleColumns = NULL) {
  # Create a function for use by resample().
  # This is called by bootstrap* and permutationTest* to create a function
  # to pass to resample, with calling sequence:
  #   resampleFun(data, ii)
  #
  # Args:
  #   data:      vector, matrix, or data frame
  #   statistic: a function or expression
  #   substituteData: substitute(data)
  #   substituteStatistic: substitute(statistic)
  #   args.stat: additional arguments to pass to statistic
  #   resampleColumns: integer or character, if supplied then only these
  #       columns are permuted.
  #
  # resampleFun's parent environment is e.g. bootstrap.
  # Get value of args.stat there.

  dimData <- dim(data)
  subscriptCall <-
    IfElse(inherits(data, "data.frame"),
           ".resampleSubscriptRows(data, ii)",
           paste0("data[ii",
                  paste(rep(",", length(dimData)), collapse = " "),
                  if(length(dimData)) " drop = FALSE", "]"))
  if(!is.null(resampleColumns)) {
    # In this case overwrite the previous subscriptCall
    if(length(dimData) != 2)
      stop("resampleColumns is only supported for matrices and data frames.")
    subscriptCall <-
      IfElse(inherits(data, "data.frame"),
             paste("{.resampleData <- data;",
                   ".resampleData[resampleColumns] <-",
                   ".resampleSubscriptRows(data[resampleColumns], ii);",
                   ".resampleData}"),
             paste("{.resampleData <- data;",
                   ".resampleData[, resampleColumns] <- ",
                   "data[ii, resampleColumns, drop = FALSE];",
                   ".resampleData}"))
  }

  # There are two cases: statistic may be a function, or an expression.
  if(is.call(substituteStatistic) && substituteStatistic[[1]] != "function") {
    # 'statistic' is an expression (function call). Basic form:
    # f <- function(data, ii) {
    #   X <- data[ii, ]  # where X the original name of the data, or 'data'
    #   eval(substituteStatistic, X, parent.frame())
    # }
    if(!is.name(substituteData)) {
      substituteData <- "data"
    } else {
      # Support cases like
      #   bootstrap(myVector, mean(data))
      # where user gives expression using data instead of myVector
      temp <- all.names(substituteStatistic)
      if(!is(data, "data.frame") &&
         is.element("data", temp) &&
         !is.element(as.character(substituteData), temp))
        substituteData <- "data"
    }
    ftext <- paste0("function(data, ii) { ",
                    substituteData, "<-", subscriptCall, "\n",
                    IfElse(is(data, "data.frame"),
                           paste0("eval(Quote(", deparse(substituteStatistic),
                                  "), ", substituteData, ")"),
                           deparse(substituteStatistic)),
                    "}")
    f <- eval(parse(text = ftext)[[1]])
    environment(f) <- new.env(parent = parent.frame())
    return(f)
  }

  # 'statistic' is a function.
  # Simplest form: f <- function(data, ii) statistic(data[ii, ])
  if(is.null(args.stat)) {
    ftext <- paste0("function(data, ii) statistic(",
                    subscriptCall, ")")
  } else {
    # Handle args.stat case, using do.call.
    ftext <- paste0("function(data, ii) { ",
                    'do.call("statistic", ',
                    "c(list(", subscriptCall, "), args.stat))}")
    # Assume that args.stat exists in parent's frame (otherwise could
    # assign it to f's frame).
  }
  if(is.name(substituteStatistic))
    ftext <- sub("statistic", substituteStatistic, ftext)
  f <- eval(parse(text = ftext)[[1]])
  environment(f) <- new.env(parent = parent.frame())
  if(!is.name(substituteStatistic))
    assign("statistic", statistic, envir = environment(f))
  if(!is.null(resampleColumns))
    assign("resampleColumns", resampleColumns, envir = environment(f))
  return(f)
}


### Some informal testing code
if(FALSE) {

  # check subscriptCall
  f <- function(dimData) paste0("data[ii",
                          paste(rep(",", length(dimData)), collapse = " "),
                          if(length(dimData)) " drop = FALSE", "]")
  f(NULL)
  f(5)
  f(5:6)
  f(5:7)

  boot <- function(data, statistic, args.stat = NULL, resampleColumns = NULL) {
    f <- .resampleMakeFunction(data, statistic, substitute(data),
                               substitute(statistic), args.stat = args.stat,
                               resampleColumns = resampleColumns)
    print(f)
    envf <- ls(environment(f))
    if(length(envf)) {
      cat("environment contains:", paste(envf, collapse = " "), "\n")
    }
    print(f(data, 1:IfElse(is.null(dim(data)), length(data), nrow(data))))
    invisible(f)
  }
  source("~/resample/R/MakeFunction.R")
  x9 <- 1:9
  xDF <- data.frame(a = 1:5, b = 2:6)

  # base case: data by name, statistic is function by name
  f <- boot(x9, mean)

  # data expression
  f <- boot(1:9, mean)

  # data expression, matrix
  f <- boot(cbind(A = 1:9, B = 2:10), colMeans)

  # args.stat
  f <- boot(x9, mean, args.stat = list(trim = .25))
  f(x9, -2) # 5.5 = trimmed mean

  # inline function
  f <- boot(x9, function(z) mean(z))

  # data frame,
  f <- boot(xDF, colMeans)


  #### statistic expression
  # data by name
  f <- boot(x9, mean(x9))

  # data as expression, refer to 'data'
  f <- boot(1:9, mean(data))

  # data frame
  f <- boot(xDF, mean(a))

  # data frame expression
  f <- boot(data.frame(a = 1:5, b = 2:6), mean(a))


  # data frame, resampleColumns integer
  f <- boot(xDF, cor, resampleColumns = 2)

  # data frame, resampleColumns character
  f <- boot(xDF, cor, resampleColumns = "b")

  # data expression, matrix, resampleColumns
  f <- boot(cbind(A = 1:9, B = 2:10), cor, resampleColumns = 2)

  # data frame, resampleColumns integer, statistic expression
  f <- boot(xDF, cor(a, b), resampleColumns = 2)


  # TODO: turn those into do.test

  source("~/resample/R/MakeFunction.R")
}
