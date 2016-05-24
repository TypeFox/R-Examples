#' Convert a package worth of RUnit tests to testthat tests
#' 
#' Converts all RUnit tests in a package to testthat tests, and writes them out 
#' to a file.
#' @param pkg Either a string containing a path to a package or a 
#' \code{devtools::package} object.
#' @param test_dir Directory inside the package containing test files.
#' @param test_file_regexp Regular expression determining which files are 
#' considered to contain tests.
#' @param test_func_regexp Regular expression determining which functions in the 
#' tests files are considered to be tests.
#' @param testthat_files Character vector of paths for the output files.  
#' Defaults to \code{stdout()} to prevent you overwriting your existing test files.
#' Use \code{runit_files} to refer to the original \code{RUnit} test files.
#' @param ... Not currently used.
#' @return A list of lists of calls.  Each call element is a \code{testthat} 
#' test. The names of the top level list correspond to the names of the input 
#' test files.  The names of the sublists correspond to the names of the 
#' \code{RUnit} test functions in that file.
#' Each list element has an environment containing the original \code{RUnit}, 
#' stored in \code{attr(, "runit_tests")}.
#' @note \code{RUnit} tests will be evaluated when they are read in, so make 
#' sure that all your tests pass before you convert them.
#' @seealso \code{\link{convert_test_file}}, \code{\link{convert_test}}
#' @export
convert_package_tests <- function(pkg, test_dir = "inst/tests", 
  test_file_regexp = "^runit.+\\.[rR]", test_func_regexp = "^test.+", 
  testthat_files = stdout(), ...)
{
  UseMethod("convert_package_tests")
}

#' @method convert_package_tests package
#' @rdname convert_package_tests 
#' @export
convert_package_tests.package <- function(pkg, test_dir = "inst/tests", 
  test_file_regexp = "^runit.+\\.[rR]", test_func_regexp = "^test.+", 
  testthat_files = stdout(), ...)
{
  convert_package_tests(
    pkg$path, 
    test_dir, 
    test_file_regexp, 
    test_func_regexp, 
    testthat_files,
    ...
  )
}

#' @method convert_package_tests character
#' @rdname convert_package_tests
#' @importFrom assertive is_empty
#' @export
convert_package_tests.character <- function(pkg, test_dir = "inst/tests", 
  test_file_regexp = "^runit.+\\.[rR]", test_func_regexp = "^test.+", 
  testthat_files = stdout(), ...)
{
  runit_files <- dir(
    file.path(pkg, test_dir), 
    pattern    = test_file_regexp, 
    full.names = TRUE
  )
  if(is_empty(runit_files))
  {
    warning("There are no test files to convert.")
    return(invisible(list()))
  }
  names(runit_files) <- runit_files
  # Weird evaluation is needed so that testthat_files can be a function of 
  # runit_files.
  # e.g., testthat_files = str_replace(runit_files, "^runit", "testthat"))
  # is useful
  testthat_files <- eval(substitute(testthat_files), envir = sys.frame(sys.nframe()))
  Map(convert_test_file, runit_files, testthat_file = testthat_files)
}

#' Convert a file worth of RUnit tests to testthat tests
#' 
#' Converts all RUnit tests in a file to testthat tests, and writes them out 
#' to another file.
#' @param runit_file A path to an file containing \code{RUnit} tests.
#' @param test_func_regexp Regular expression determining which functions in the 
#' tests files are considered to be tests.
#' @param testthat_file String of path for the output file.  
#' Defaults to \code{stdout()} to prevent you overwriting your existing test 
#' files.  See note.
#' @return A list of calls to \code{test_that} is invisibly returned.  This list
#' is also written to an output connection (defaulting to stdout).
#' An environment containing the original \code{RUnit} is stored in 
#' \code{attr(, "runit_tests")}.
#' @note \code{RUnit} tests will be evaluated when they are read in, so make 
#' sure that all your tests pass before you convert them.
#' 
#' The \code{testthat_file} argument can be based upon the \code{runit_file}
#' argument.  For example, if you use the traditional \code{RUnit} test file 
#' naming strategy, something like \code{sub("^runit", "testthat", runit_file)}
#' may be appropriate.
#' @importFrom assertive is_empty
#' @importFrom assertive is_stdout
#' @seealso \code{\link{convert_package_tests}}, \code{\link{convert_test}}
#' @examples
#' tmp <- tempfile("test-convert_test_file")
#' writeLines(
#'   "test_truth <- function()
#' {
#'   x <- all(runif(10) > 0)
#'   checkTrue(x)
#' }
#' test_equality <- function()
#' {
#'   x <- sqrt(1:5)
#'   expected <- c(1, 4, 9, 16, 25)
#'   checkEquals(expected, x ^ 4)
#' }
#' test_error <- function()
#' {
#'   checkException('1' + '2')
#' }",
#'   tmp
#' )
#' convert_test_file(tmp)
#' unlink(tmp)
#' @export
convert_test_file <- function(runit_file, test_func_regexp = "^test.+", 
  testthat_file = stdout())
{
  e <- new.env()
  sys.source(runit_file, envir = e)
  test_fn_names <- ls(e, pattern = test_func_regexp)
  if(is_empty(test_fn_names))
  {
    warning(
      "There are no test functions in the file ", 
      get_filename(runit_file), 
      " to convert."
    )
  }
  names(test_fn_names) <- test_fn_names
  new_tests <- lapply(
    test_fn_names,
    function(test_fn_name)
    {
      convert_test(
        get(test_fn_name, envir = e), 
        test_fn_name
      )
    }
  )
  attr(new_tests, "runit_tests") <- e
  output <- unlist(
    lapply(
      new_tests,
      function(new_test)
      {
        c(deparse(new_test), "")
      }
    )  
  )
  if(is_stdout(testthat_file) || 
     !file.exists(as.character(testthat_file)) || 
     get_yes_no_response("The file ", get_filename(testthat_file), "already exists. Overwrite (yes/no)? > "))
  {
    message("Writing tests to ", get_filename(testthat_file))
    cat(output, file = testthat_file, sep = "\n")
  }  
  invisible(new_tests)
}

#' Convert an RUnit test to a testthat test
#' 
#' Converts a single RUnit test to a testthat test.
#' @param runit_test_fn An \code{RUnit} test function.
#' @param test_description A string describing the test, defaulting to
#' the name of the \code{RUnit} test function.
#' @return A call to \code{test_that}, containing a \code{testhat} test 
#' equivalent to the input \code{RUnit} test.
#' @seealso \code{\link{convert_package_tests}}, \code{\link{convert_test_file}}
#' @examples
#' test_truth <- function()
#' {
#'   x <- all(runif(10) > 0)
#'   checkTrue(x)
#' }
#' convert_test(test_truth)
#' test_equality <- function()
#' {
#'   x <- sqrt(1:5)
#'   expected <- c(1, 4, 9, 16, 25)
#'   checkEquals(expected, x ^ 4)
#' }
#' convert_test(test_equality)
#' test_error <- function()
#' {
#'   checkException("1" + "2")
#' }
#' convert_test(test_error)
#' @export
convert_test <- function(runit_test_fn, 
  test_description = assertive::get_name_in_parent(runit_test_fn))
{
  new_test_body <- convert_lines(body(runit_test_fn))
  call(
    "test_that",
    test_description,
    new_test_body
  )  
}

#' Convert lines of RUnit code to testthat code
#' 
#' Converts lines of \code{RUnit} code to the equivalent \code{testthat} code.
#' @param lines An object of class \code{call}.
#' @return An object of class \code{call}.
convert_lines <- function(lines)
{
  lines <- as.list(lines)
  for(i in seq_along(lines))
  {
    lines[[i]] <- convert_line(lines[[i]])    
  } 
  as.call(lines)
}

#' Convert lines of RUnit code to testthat code
#' 
#' Converts lines of \code{RUnit} code to the equivalent \code{testthat} code.
#' @param x An R language object.
#' @return An object of class \code{call}.
convert_line <- function(x)
{
  current_fn_name <- deparse(as.list(x)[[1]])
  switch(
    current_fn_name,
    checkEquals        = convert_checkEquals(x),
    checkEqualsNumeric = convert_checkEqualsNumeric(x), 
    checkException     = convert_checkException(x),
    checkIdentical     = convert_checkIdentical(x),
    checkTrue          = convert_checkTrue(x),
    "for"              = convert_for(x),
    "while"            = convert_while(x),
    "repeat"           = convert_repeat(x),
    "if"               = convert_if(x),
    "switch"           = convert_switch(x),
    "{"                = convert_brace(x),
    x
  )
}

#' Convert a checkEquals call.
#' 
#' Converts a \code{checkEquals} call.
#' 
#' @param x A call to the \code{\link[RUnit]{checkEquals}} function.
#' @return x A call to the \code{\link[testthat]{expect_equal}} function.
convert_checkEquals <- function(x)
{
  the_line <- as.list(match.call(RUnit::checkEquals, x))
  fixup(
    call(
      "expect_equal", 
      the_line$current, 
      the_line$target,
      info = the_line$msg
    )
  )
}

#' Convert a checkEqualsNumeric call.
#' 
#' Converts a \code{checkEqualsNumeric} call.
#' 
#' @param x A call to the \code{\link[RUnit]{checkEqualsNumeric}} function.
#' @return x A call to the \code{\link[testthat]{expect_equal}} function.
convert_checkEqualsNumeric <- function(x)
{
  the_line <- as.list(match.call(RUnit::checkEqualsNumeric, x))
  fixup(
    call(
      "expect_equal", 
      the_line$current, 
      the_line$target,
      info = the_line$msg
    )
  )
}

#' Convert a checkException call.
#' 
#' Converts a \code{checkException} call.
#' 
#' @param x A call to the \code{\link[RUnit]{checkException}} function.
#' @return x A call to the \code{\link[testthat]{expect_error}} function.
convert_checkException <- function(x)  
{
  the_line <- as.list(match.call(RUnit::checkException, x))
  fixup(call(
    "expect_error", 
    the_line$expr,
    info = the_line$msg
  ))
}  

#' Convert a checkIdentical call.
#' 
#' Converts a \code{checkIdentical} call.
#' 
#' @param x A call to the \code{\link[RUnit]{checkIdentical}} function.
#' @return x A call to the \code{\link[testthat]{expect_identical}} function.
convert_checkIdentical <- function(x)
{
  the_line <- as.list(match.call(RUnit::checkIdentical, x))
  fixup(call(
    "expect_identical", 
    the_line$current, 
    the_line$target,
    info = the_line$msg
  ))
}  

#' Convert a checkTrue call.
#' 
#' Converts a \code{checkTrue} call.
#' 
#' @param x A call to the \code{\link[RUnit]{checkTrue}} function.
#' @return x A call to the \code{\link[testthat]{expect_true}} or 
#' \code{\link[testthat]{expect_false}} function. (The later occurs if the input
#' to \code{checkTrue} was negated.)
convert_checkTrue <- function(x)
{
  the_line <- as.list(match.call(RUnit::checkTrue, x))
  the_object <- as.list(the_line$expr)      
  if(identical(deparse(the_object[[1]]), "!") && identical(match.fun(the_object[[1]]), `!`))
  {
    fixup(
      call(
        "expect_false", 
        the_object[[2]],
        info = the_line$msg
      )
    )
  } else
  {
    fixup(
      call(
        "expect_true",
        the_line$expr,
        info = the_line$msg
      )
    )
  }      
}

#' Convert a for block.
#' 
#' Converts a \code{for} call.
#' 
#' @param x A call to the \code{\link[base]{for}} function.
#' @return x An object of class \code{call}.
convert_for <- function(x)
{
  the_loop <- as.list(x)
  call("for", the_loop[[2]], the_loop[[3]], convert_line(the_loop[[4]]))
}

#' Convert a while block.
#' 
#' Converts a \code{while} call.
#' 
#' @param x A call to the \code{\link[base]{while}} function.
#' @return x An object of class \code{call}.
convert_while <- function(x)
{
  the_loop <- as.list(x)
  call("while", the_loop[[2]], convert_lines(the_loop[[3]]))
}

#' Convert a repeat block.
#' 
#' Converts a \code{repeat} call.
#' 
#' @param x A call to the \code{\link[base]{repeat}} function.
#' @return x An object of class \code{call}.
convert_repeat <- function(x)
{
  the_loop <- as.list(x)
  call("repeat", convert_lines(the_loop[[2]]))
}

#' Convert an if block.
#' 
#' Converts an \code{if} call.
#' 
#' @param x A call to the \code{\link[base]{if}} function.
#' @return x An object of class \code{call}.
convert_if <- function(x)
{
  the_flow <- as.list(x)
  if(length(the_flow) == 3)
  {
    call("if", the_flow[[2]], convert_lines(the_flow[[3]]))
  } else
  {
    # there is an else clause
    call("if", the_flow[[2]], convert_lines(the_flow[[3]]), convert_lines(the_flow[[4]])) 
  }
}

#' Convert a switch block.
#' 
#' Converts a \code{switch} call. 
#' 
#' @param x A call to the \code{\link[base]{switch}} function.
#' @return x the input, with a warning.
convert_switch <- function(x)
{
  the_flow <- as.list(x)
  converted <- lapply(the_flow[-1:-2], convert_line)
  as.call(c(as.name("switch"), the_flow[[2]], converted))
}

#' Convert a code block.
#' 
#' Converts a \code{switch} call.  
#' 
#' @param x A call to the brace function.
#' @return x the input, with a warning.
convert_brace <- function(x)
{
  the_block <- as.list(x)
  if(length(the_block) == 1)
  {
    return(x)
  } else 
  {
    converted <- lapply(the_block[-1], convert_line)
    as.call(c(as.name("{"), converted))
  }
}

#' Get a quoted filename
#' 
#' Gets a single quoted filename from a connection or path.
#' @param x A \code{connection} object or a character vector of file paths.
#' @return A character vector of paths to files.
get_filename <- function(x)
{
  sQuote(
    if(inherits(x, "connection"))
    {
      summary(x)$description
    } else
    {
      as.character(x)
    }
  )
}
