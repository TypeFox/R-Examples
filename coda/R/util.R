"read.yesno" <-
function (string, default=TRUE)
{
  wrd <- ifelse(default, " (Y/n)?\n:", " (y/N)?\n:")
  cat("\n", string, wrd, sep = "")
  ans <- readline()
  val <- if (default) 
    pmatch(ans, c("no","NO"), nomatch=0) == 0
  else
    pmatch(ans, c("yes","YES"), nomatch=0) != 0
  return(val)
}

"change.tfoption" <-
function (string, option) 
{
  current.value <- coda.options(option)
  if (!is.logical(current.value)) 
    stop("Invalid option: must take logical values")
  new.value <- read.yesno(string, current.value)
  if (new.value != current.value) {
    arg <- list(new.value)
    names(arg) <- option
    coda.options(arg)
  }
  invisible()
}

"coda.options" <-
function (...) 
{
  ## Set and display coda options
  single <- FALSE
  copt <- if (exists("options", envir=coda.env, inherits=FALSE)) {
      get("options", envir=coda.env, inherits=FALSE)
  }
  else {
    .Coda.Options.Default
  }
  if (nargs() == 0) {
    return(copt)
  }
  else {
    args <- list(...)
    if (length(args) == 1) {
      if (is.list(args[[1]])) 
        args <- args[[1]]
      else if (is.null(names(args))) 
        single <- TRUE
    }
  }
  if (is.null(names(args))) {
    ## Display options
    args <- unlist(args)
    value <- vector("list", length(args))
    names(value) <- args
    for (v in args) if (any(v == names(copt))) 
      value[v] <- copt[v]
    if (single) 
      return(value[[1]])
    else return(value)
  }
  else {
    ## Set options
    oldvalue <- vector("list", length(args))
    names(oldvalue) <- names(args)
    if (any(names(args) == "default") && args$default == 
        TRUE) 
      copt <- .Coda.Options.Default
    for (v in names(args)) if (any(v == names(copt))) {
      oldvalue[v] <- copt[v]
      if (is.null(args[[v]])) 
        copt[v] <- list(NULL)
      else if (mode(copt[[v]]) == mode(args[[v]])) 
        copt[v] <- args[v]
    }
    assign("options", copt, envir=coda.env)
    invisible(oldvalue)
  }
}

"multi.menu" <- function (choices, title, header, allow.zero = TRUE) 
{
  ## Select more than one value from a menu 
  ## 

  if (!missing(title)) 
    cat(title, "\n\n")
  mat <- matrix(c(1:length(choices), choices), ncol = 2)
  if (!missing(header)) {
    if (length(header) == 2) 
      mat <- rbind(header, mat)
    else stop("header is wrong length")
  }
  cat(paste(format(mat[, 1]), format(mat[, 2])), sep = "\n")
  repeat {
    cat("\nEnter relevant number(s), separated by commas", 
        "Ranges such as 3:7 may be specified)", sep = "\n")
    if (allow.zero) 
      cat("(Enter 0 for none)\n")
    ans <- scan(what = character(), sep = ",", strip.white = TRUE, 
                nlines = 1, quiet = TRUE)
    if (length(ans) > 0) {
      out <- numeric(0)
      for (i in 1:length(ans)) {
        nc <- nchar(ans[i])
        wrd <- substring(ans[i], 1:nc, 1:nc)
        colons <- wrd == ":"
        err <- any(is.na(as.numeric(wrd[!colons]))) | 
        sum(colons) > 1 | colons[1] | colons[nc]
        if (err) {
          cat("Error: you have specified a non-numeric value!\n")
          break
        }
        else {
          out <- c(out, eval(parse(text = ans[i])))
          if (min(out) < ifelse(allow.zero, 0, 1) | max(out) > 
              length(choices) | (any(out == 0) & length(out) > 
                                 1)) {
            err <- TRUE
            cat("Error: you have specified variable number(s) out of range!\n")
            break
          }
        }
      }
      if (!err) 
        break
    }
  }
  return(out)
}

"read.and.check" <-
  function (message = "", what = numeric(), lower, upper, answer.in, default) 
{
  ## Read data from the command line and check that it satisfies 
  ## certain conditions.  The function will loop until it gets 
  ## and answer satisfying the conditions. This entails extensive 
  ## checking of the conditions to  make sure they are consistent 
  ## so we don't end up in an infinite loop. 

  have.lower <- !missing(lower)
  have.upper <- !missing(upper)
  have.ans.in <- !missing(answer.in)
  have.default <- !missing(default)
  if (have.lower | have.upper) {
    if (!is.numeric(what)) 
      stop("Can't have upper or lower limits with non numeric input")
    if (have.lower && !is.numeric(lower)) 
      stop("lower limit not numeric")
    if (have.upper && !is.numeric(upper)) 
      stop("upper limit not numeric")
    if ((have.upper & have.lower) && upper < lower) 
      stop("lower limit greater than upper limit")
  }
  if (have.ans.in) {
    if (mode(answer.in) != mode(what)) 
      stop("inconsistent values of what and answer.in")
    if (have.lower) 
      answer.in <- answer.in[answer.in >= lower]
    if (have.upper) 
      answer.in <- answer.in[answer.in <= upper]
    if (length(answer.in) == 0) 
      stop("No possible response matches conditions")
  }
  if (have.default) {
    if (mode(default) != mode(what)) 
      stop("inconsistent values of what and default")
    if (have.lower && default < lower) 
      stop("default value below lower limit")
    if (have.upper && default > upper) 
      stop("default value above upper limit")
    if (have.ans.in && !any(answer.in == default)) 
      stop("default value does not satisfy conditions")
  }
  err <- TRUE
  while (err) {
    if (nchar(message) > 0) {
      cat("\n", message, "\n", sep = "")
      if (have.default) 
        cat("(Default = ", default, ")\n", sep = "")
    }
    repeat {
      cat("1:")
      ans <- readline()
      if (length(ans) == 1 && nchar(ans) > 0) 
        break
      else if (have.default) {
        ans <- default
        break
      }
    }
    if (is.numeric(what)) {
      err1 <- TRUE
      ans <- as.numeric(ans)
      message <- "You must enter a number"
      if (is.na(ans)) 
        NULL
      else if ((have.lower & have.upper) && (ans < lower | 
                                             ans > upper)) 
        message <- paste(message, "between", lower, "and", 
                         upper)
      else if (have.lower && ans < lower) 
        message <- paste(message, ">=", lower)
      else if (have.upper && ans > upper) 
        message <- paste(message, "<=", upper)
      else err1 <- FALSE
    }
    else err1 <- FALSE
    if (have.ans.in) {
      if (!is.na(ans) && !any(ans == answer.in)) {
        message <- paste("You must enter one of the following:", 
                         paste(answer.in, collapse = ","))
        err2 <- TRUE
      }
      else err2 <- FALSE
    }
    else err2 <- FALSE
    err <- err1 | err2
  }
  return(ans)
}

".Coda.Options.Default" <-
  list(trace = TRUE,
       densplot = TRUE,
       lowess = FALSE, 
       combine.plots = TRUE,
       bandwidth = function (x) 
       {
         x <- x[!is.na(x)]
         1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2
       },
       digits = 3,
       quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
       frac1 = 0.1,
       frac2 = 0.5,
       q = 0.025,
       r = 0.005, 
       s = 0.95,
       combine.stats = FALSE,
       combine.corr = FALSE,
       halfwidth = 0.1,
       user.layout = FALSE,
       gr.bin = 10,
       geweke.nbin = 20,
       gr.max = 50
       )

coda.env <- new.env()











