#' bisectr package
#'
#' This package is used for creating test scripts to find bad commits with
#' git bisect.
#' For example test scripts, see \url{https://github.com/wch/bisectr}.
#'
#' @name bisectr
#' @docType package
#' @aliases bisectr package-bisectr
NULL


#' Run a test function for git bisect testing.
#'
#' If the function \code{fun} returns \code{"good"} or \code{TRUE},
#' quit and return a code to mark this commit as good.
#' If the function returns \code{"bad"} or \code{FALSE},
#' quit and return a code to mark this commit as bad.
#' If the function returns \code{"skip"} or \code{NA},
#' quit and return a code to mark this commit as skip.
#' If the function returns \code{"ignore"} or \code{NULL}, do nothing.
#'
#' It is also important to set \code{on_error}. This tells it what to
#' do when the test function throws an error. The default behavior is to
#' mark this commit as skip. However, in some cases, it makes
#' sense to mark this commit as bad if an error is thrown.
#'
#' @seealso \code{\link{bisect_load_all}}
#' @seealso \code{\link{bisect_install}}
#' @seealso \code{\link{bisect_source}}
#' @seealso \code{\link{bisect_return_interactive}}
#' 
#' @param fun      The test function
#' @param on_error What to do if running \code{fun} throws an error 
#'                  (default is to mark this commit as skip)
#' @param msg  A message to print to the console when running the test
#' @export
bisect_runtest <- function(fun, on_error = "skip", msg = "Running test...") {

  # Check that fun is a function -- easy to accidentally pass myfun()
  # instead of myfun.
  if (!is.function(fun)) {
    stop("'fun' is not a function. Make sure to pass 'myfunction' and not 'myfunction()'")
  }

  message(msg)

  error_fun <- function(e) {
    message(e)
    message("\nError encountered in test.")
    return(on_error)
  }

  status <- tryCatch(fun(), error = error_fun)
  
  # The identical() bit is necessary so that NULL and NA comparisons work
  if (is.null(status) || identical(tolower(status), "ignore")) {
    # Return NULL, but don't print
    invisible(NULL)
  } else if (is.na(status) || identical(tolower(status), "skip")) {
    mark_commit_skip()
  } else if (status == TRUE || status == "good") {
    mark_commit_good()
  } else if (status == FALSE || status == "bad") {
    mark_commit_bad()
  }
}


#' Like \code{source}, but for bisect tests.
#'
#' If the file fails to load, the default is mark this commit as skip.
#'
#' @seealso \code{\link{source}}
#' @seealso \code{\link{bisect_load_all}}
#' @seealso \code{\link{bisect_install}}
#' @seealso \code{\link{bisect_runtest}}
#' @seealso \code{\link{bisect_return_interactive}}
#'
#' @param file     The file to load
#' @param ...      Other arguments to pass to \code{\link{source}}
#' @param on_error What to do if loading throws an error (default is to mark this
#'  commit as "skip")
#' @export
#' @importFrom devtools load_all
bisect_source <- function(file, ..., on_error = "skip") {
  bisect_runtest(function() {
      source(file, ...)
      return("good")
    },
    on_error = on_error,
    msg = paste("Sourcing file ", file))
}


#' Like \code{load_all}, but for bisect tests.
#'
#' If the package fails to load, the default is to mark this commit as skip.
#'
#' @seealso \code{\link{bisect_source}}
#' @seealso \code{\link{bisect_install}}
#' @seealso \code{\link{bisect_runtest}}
#' @seealso \code{\link{bisect_return_interactive}}
#'
#' @param pkgdir   The directory to load from
#' @param on_error What to do if loading throws an error (default is to mark this
#'  commit as "skip")
#' @export
#' @importFrom devtools load_all
bisect_load_all <- function(pkgdir = ".", on_error = "skip") {
  bisect_runtest(function() {
      load_all(pkgdir, reset = TRUE)
    }, 
    on_error = on_error,
    msg = paste("Loading package in directory", pkgdir))
}


#' Install a package from source, for bisect tests.
#'
#' If the installation fails, the default behavior is to mark this commit
#' as skip.
#'
#' This function is usually used together with \code{bisect_require}.
#'
#' @seealso \code{\link{bisect_require}}
#' @seealso \code{\link{bisect_load_all}}
#' @seealso \code{\link{bisect_source}}
#' @seealso \code{\link{bisect_runtest}}
#' @seealso \code{\link{bisect_return_interactive}}
#'
#' @param pkgdir  The directory to load from
#' @param on_fail What to do if installation fails (default is to mark this
#'  commit as "skip")
#' @export
#' @importFrom devtools dev_mode
#' @importFrom devtools install
bisect_install <- function(pkgdir = ".", on_fail = "skip") {
  tempPkgdir <- normalizePath(file.path(tempdir(), "bisect-pkgs"),
                              winslash = "/", mustWork = FALSE)
  dev_mode(TRUE, path = tempPkgdir)
  message("Temp package installation directory: ", tempPkgdir)

  # install() returns TRUE on success; in this case, we'll give a "ignore" code
  #   so that the test script will continue.
  # When install() fails, it throws an error, in which case we'll pass along
  #  the on_fail code.
  bisect_runtest(function() {
      install(pkgdir)
      return("ignore")
    },
    on_error = on_fail,
    msg = paste("Installing package in directory", pkgdir)
  )
}


#' Load a package like \code{require()}, for bisect tests.
#'
#' If the package fails to load, the default behavior is to mark this commit
#' as skip.
#'
#' This function is usually used together with \code{bisect_install}.
#'
#' @seealso \code{\link{bisect_install}}
#' @seealso \code{\link{bisect_load_all}}
#' @seealso \code{\link{bisect_source}}
#' @seealso \code{\link{bisect_runtest}}
#' @seealso \code{\link{bisect_return_interactive}}
#'
#' @param package Name of package
#' @param on_fail What to do if loading fails (default "skip")
#' @export
bisect_require <- function(package, on_fail = "skip") {
  
  package <- as.character(substitute(package))

  # With require(), success returns TRUE and failure returns FALSE
  # but we need to pass different return values to bisect_runtest().
  # If success loading, do nothing ("ignore"); if failure, return on_fail
  bisect_runtest(function() {
      if (require(package, character.only = TRUE))
        return("ignore")
      else
        return(on_fail)
    },
    msg = paste("Loading package", package)
  )
}


#' Prompt the user for an interactive good/bad/skip response and return
#' the appropriate value (to be passed to \code{bisect_runtest}).
#'
#' @seealso \code{\link{bisect_runtest}}
#' @seealso \code{\link{bisect_load_all}}
#' @seealso \code{\link{bisect_install}}
#' @seealso \code{\link{bisect_source}}
#'
#' @export
bisect_return_interactive <- function() {
  while (1) {
    message("Mark this commit [g]ood, [b]ad, or [s]kip? ", appendLF = FALSE)
    
    # Need to use "stdin" to get user input in a script -- stdin() doesn't work
    response <- scan("stdin", what = character(), n = 1, quiet = TRUE) 
    
    if (identical(tolower(response), "g")) {
      return("good")
    } else if (identical(tolower(response), "b")) {
      return("bad")
    } else if (identical(tolower(response), "s")) {
      return("skip")
    } else {
      message(paste("Unknown response:", response))
    }
  }
}


# ===========================================================
# Functions to quit with return code for marking commits

mark_commit_good <- function() {
  message("Returning code: good (0)\n")
  quit(status = 0)
}

mark_commit_bad <- function() {
  message("Returning code: bad (1)\n")
  quit(status = 1)
}

mark_commit_skip <- function() {
  message("Returning code: skip (125)\n")
  quit(status = 125)
}
