#' @rdname r_can_find_tools
#' @export
assert_r_can_find_tools <- function(tools, severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_can_find_tools, tools = tools, severity = severity)
}

#' @rdname r_can_find_tools
#' @export
assert_r_can_compile_code <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_can_compile_code, severity = severity)
}

#' @rdname r_can_find_tools
#' @export
assert_r_can_build_translations <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_can_build_translations, severity = severity)
}

#' @rdname r_can_find_tools
#' @export
assert_r_can_find_java <- function(java_type = c("same_as_r", "any", "64bit", "32bit"), severity = getOption("assertive.severity", "stop"))
{
  java_type <- match.arg(java_type)
  assert_engine(r_can_find_java, java_type = java_type, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_jpeg_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_jpeg_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_png_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_png_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_tiff_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_tiff_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_tcltk_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_tcltk_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_x11_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_x11_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_aqua_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_aqua_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_http_ftp_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_http_ftp_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_sockets_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_sockets_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_libxml_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_libxml_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_fifo_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_fifo_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_cledit_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_cledit_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_iconv_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_iconv_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_nls_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_nls_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_profmem_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_profmem_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_cairo_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_cairo_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_icu_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_icu_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_long_double_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_long_double_capability, severity = severity)
}

#' @rdname r_has_jpeg_capability
#' @export
assert_r_has_libcurl_capability <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(r_has_libcurl_capability, severity = severity)
}
