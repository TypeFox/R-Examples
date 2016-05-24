#' @include imports.R

#' @rdname is_connection
#' @export
assert_is_bzfile_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_bzfile_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_fifo_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_fifo_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_file_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_file_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  ) 
}

#' @rdname is_connection
#' @export
assert_is_gzfile_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_gzfile_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_incomplete_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_incomplete_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_open_connection <- function(x, rw = "", 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_open_connection, 
    x, 
    rw = rw, 
    .xname = get_name_in_parent(x),
    severity = severity
  )   
}

#' @rdname is_connection
#' @export
assert_is_pipe_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_pipe_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_readable_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_readable_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_socket_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_socket_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_stderr <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_stderr, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_stdin <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_stdin, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_stdout <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_stdout, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_terminal_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_terminal_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_text_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_text_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )   
}

#' @rdname is_connection
#' @export
assert_is_unz_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_unz_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_url_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_url_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_connection
#' @export
assert_is_writable_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_writable_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  ) 
}

#' @rdname is_connection
#' @export
assert_is_xzfile_connection <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_xzfile_connection, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
