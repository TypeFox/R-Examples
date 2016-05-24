#' Query construction
#'
#' @export
#' @param .data Result of a call to \code{api}
#' @param ...	Comma separated list of unquoted variable names. These are combined into
#' a list and passed to whatever http method is used downstream
#' @param .dots	Used to work around non-standard evaluation
#' @param body_value one of the following:
#' \itemize{
#'  \item FALSE: No body
#'  \item NULL: An empty body
#'  \item "": A length 0 body
#'  \item upload_file("path/"): The contents of a file. The mime type will be guessed
#'  from the extension, or can be supplied explicitly as the second argument
#'  to upload_file()
#'  \item A character or raw vector: sent as is in body. Use content_type to tell the
#'  server what sort of data you are sending.
#' }
#' @family dsl
#' @examples \dontrun{
#' ## NSE
#' dd <- api("http://httpbin.org/post")
#' dd %>% api_body(body_value = NULL) %>% http("POST")
#' dd %>% api_body(body_value = "") %>% http("POST")
#'
#' ## other named parameters are passed as form values
#' dd %>% api_body(x = hello) %>% http("POST")
#'
#' # upload a file
#' file <- "~/some_test.txt"
#' cat("hello, world", file = file)
#' dd %>% api_body(x = upload_file("~/some_test.txt")) %>% http("POST")
#'
#' # A named list
#' dd %>% api_body(x = hello, y = stuff) %>% http("POST")
#'
#' ## SE
#' dd %>% api_body_(x = "hello", y = "stuff") %>% http("POST")
#' }
api_body <- function(.data, ..., body_value = NULL){
  api_body_(.data, .dots = lazyeval::lazy_dots(...), body_value = body_value)
}

#' @export
#' @rdname api_body
api_body_ <- function(.data, ..., .dots, body_value = NULL){
  ## FIXME - need to toggle on POST by default here when body passed
  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  if (is.null(body_value)) {
    dots <- lazyeval::all_dots(.dots, ...)
    args <- sapply(dots, "[[", "expr")
    modifyList(.data, list(body = as.list(args)))
  } else {
    modifyList(.data, list(body = body_value))
  }
}
