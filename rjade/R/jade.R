#' Render Jade Template
#'
#' Jade is a high performance template engine heavily influenced by Haml.
#'
#' Converting a template to HTML text involves two steps. The first step compiles
#' the template with some formatting options into a closure. The binding for this
#' is implemented in \code{jade_compile}. The second step calls the closure with
#' optionally some local variables to render the output to HTML.
#'
#' The \code{jade_render} function is a convenience wrapper that does both steps at
#' once. This is slightly faster if you only need to render your template once.
#'
#' @export
#' @rdname jade
#' @aliases jade
#' @references Jade documentation: \url{http://jade-lang.com}
#' @param text string with jade template.
#' @param ... options passed to the compiler, see \url{http://jade-lang.com/api}.
#' @param locals local variables used in the template.
#' @examples # Example from http://jade-lang.com
#' text <- readLines(system.file("examples/test.jade", package = "rjade"))
#'
#' # Compile and render seperately
#' tpl <- jade_compile(text, pretty = TRUE)
#' tpl()
#' tpl(youAreUsingJade = TRUE)
#'
#' # Slightly faster for one-time rendering
#' jade_render(text, pretty = TRUE)
#' jade_render(text, pretty = TRUE, locals = list(youAreUsingJade = TRUE))
jade_compile <- function(text, ...){
  text <- paste(text, collapse = "\n")
  opts <- list(...)
  ct <- new_jade()
  ct$call("jade_compile", text, opts)
  fn <- function(...){
    locals <- list(...)
    html <- ct$call("fn", locals)
    class(html) <- "jade_html"
    html
  }
  class(fn) <- c("jade", class(fn))
  fn
}

#' @export
#' @rdname jade
jade_render <- function(text, ..., locals = list()){
  text <- paste(text, collapse = "\n")
  opts <- list(...)
  html <- ct$call("jade.render", text, c(opts, locals))
  class(html) <- "jade_html"
  html
}

#' @export
print.jade <- function(x, ...){
  cat("Compiled Jade Template.\nCall me as a function to render.\n")
}

#' @export
print.jade_html <- function(x, ...){
  cat(x, "\n")
}

# Create new context with Jade
new_jade <- function(){

  # Create context
  myct <- new_context();

  # Temporary disabled due to UTF8 problems in windows. Fixed in V8 0.6
  # ct$source(system.file("js/jade.min.js", package = pkgname))

  src <- readLines(system.file("js/jade.min.js", package = "rjade"), encoding = "UTF-8", warn = FALSE)
  myct$eval(src)
  myct$source(system.file("js/bindings.js", package = "rjade"))
  return(myct)
}
