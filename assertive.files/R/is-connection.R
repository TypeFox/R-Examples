#' @rdname is_connection
#' @export
is_bzfile_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "bzfile", .xname)
}

#' Is the input a connection?
#'
#' Various checks to see if the input is a (particular type of/open/incomplete) 
#' connection.
#' @param x Input to check.
#' @param rw Read-write status of connection.  Passed to \code{isOpen}.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_connection} checks for objects of class "connection".
#' \code{is_open_connection} and \code{is_incomplete_connection} wrap 
#' \code{isOpen} and \code{isIncomplete} respectively, providing more 
#' information on failure.
#' \code{is_readable_connection} and \code{is_writable_connection} tell you
#' whether the connection is readable from or writable to.
#' \code{is_bzfile_connection}, \code{is_fifo_connection}, 
#' \code{is_file_connection}, \code{is_pipe_connection}, 
#' \code{is_socket_connection}, \code{is_stderr}, \code{is_stdin}, 
#' \code{is_stdout}, \code{is_text_connection}, \code{is_unz_connection},
#' \code{is_url_connection} and \code{is_xzfile_connection} give more
#' specific tests on the type of connection.
#' The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @note \code{is_incomplete_connection} will return false for closed 
#' connections, regardless of whether or not the connection ends with a newline 
#' character.
#' (\code{isIncomplete} throws an error for closed connections.)
#' @seealso \code{\link[base]{isOpen}}.
#' @examples
#' assert_is_terminal_connection(stdin())
#' assert_is_readable_connection(stdin())
#' assert_is_open_connection(stdin())
#' assert_is_stdin(stdin())
#' # Next line is usually true but, e.g., devtools::run_examples overrides it
#' assertive.base::dont_stop(assert_is_terminal_connection(stdout()))
#' assert_is_writable_connection(stdout())
#' assert_is_open_connection(stdout())
#' assert_is_stdout(stdout())
#' assert_is_terminal_connection(stderr())
#' assert_is_writable_connection(stderr())
#' assert_is_open_connection(stderr())
#' assert_is_stderr(stderr())
#' tcon <- textConnection("txt", "w", local = TRUE)
#' assert_is_text_connection(tcon)
#' assert_is_open_connection(tcon)
#' cat("this has no final newline character", file = tcon)
#' assert_is_incomplete_connection(tcon)
#' close(tcon)
#' # These examples should fail.
#' assertive.base::dont_stop({
#'   assert_is_connection("not a connection")
#'   assert_is_readable_connection(stdout())
#'   assert_is_writable_connection(stdin())
#' })
#' \dontrun{
#' fcon <- file()
#' close(fcon)
#' assert_is_open_connection(fcon)
#' }
#' @export
is_connection <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "connection", .xname)
}

#' @rdname is_connection
#' @export
is_fifo_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "fifo", .xname)
}

#' @rdname is_connection
#' @export
is_file_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "file", .xname)
}

#' @rdname is_connection
#' @export
is_gzfile_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "gzfile", .xname)
}

#' @rdname is_connection
#' @export
is_incomplete_connection <- function(x, .xname = get_name_in_parent(x))
{  
  if(!(ok <- is_open_connection(x))) 
  {
    return(ok)
  }
  if(!isIncomplete(x))
  {
    return(false("The connection %s is complete.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @importFrom assertive.base use_first
#' @export
is_open_connection <- function(x, rw = "", .xname = get_name_in_parent(x))
{
  rw <- use_first(rw)
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  conn_is_open <- try(isOpen(x, rw))
  if(inherits(conn_is_open, "try-error"))
  {
    return(false("The connection %s is not open.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @export
is_pipe_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  # Can't use is2(x, "pipe", .xname).  The class is OS dependent.
  summary_of_x <- summary(x)
  if(summary_of_x$class != "pipe")
  {
    return(false("%s is not a pipe", .xname))
  }
  TRUE
}
#' @rdname is_connection
#' @export
is_readable_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  summary_of_x <- summary(x)
  if(summary_of_x$`can read` != "yes")
  {
    return(false("The connection %s is not readable.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @export
is_socket_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "sockconn", .xname)
}

#' @rdname is_connection
#' @export
is_stderr <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  if(!identical(x, stderr()))
  {
    return(false("The connection %s is not stderr.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @export
is_stdin <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  if(!identical(x, stdin()))
  {
    return(false("The connection %s is not stdin.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @export
is_stdout <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  # Can't just test summary(x)$description because, e.g., devtools::run_examples
  # overrides this
  if(!identical(x, stdout()))
  {
    return(false("The connection %s is not stdout.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @export
is_terminal_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "terminal", .xname)
}

#' @rdname is_connection
#' @export
is_text_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "textConnection", .xname)
}

#' @rdname is_connection
#' @export
is_unz_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "unz", .xname)
}

#' @rdname is_connection
#' @export
is_url_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "url", .xname)
}

#' @rdname is_connection
#' @export
is_writable_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  summary_of_x <- summary(x)
  if(summary_of_x$`can write` != "yes")
  {
    return(false("The connection %s is not writable.", .xname))
  }
  TRUE
}

#' @rdname is_connection
#' @export
is_xzfile_connection <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_connection(x))) 
  {
    return(ok)
  }
  is2(x, "xzfile", .xname)
}
