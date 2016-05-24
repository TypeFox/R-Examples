#' Bootstrap 2 components for use with Shiny
#'
#' @name shinybootstrap2-package
#' @aliases shinybootstrap2
#' @docType package
#' @import htmltools jsonlite
NULL


#' Run Shiny UI code with Bootstrap 2 elements.
#'
#' This function takes an expression containing calls to functions in the
#' shinybootstrap2 package, and evaluates it in an environment where these
#' functions will be found, even when shinybootstrap2 is not attached.
#'
#' Shiny version 0.11 and above uses Bootstrap 3 instead of Bootstrap 2. The
#' purpose of the shinybootstrap2 package is to provide backward compatibility
#' when needed. Almost all of the functions in shinybootstrap2 have the same
#' name as functions in shiny, but they generate HTML that works with Bootstrap
#' 2 instead of 3.
#'
#' This function should almost always be called using
#' \code{shinybootstrap2::withBootstrap2()}, without attaching the package. In
#' other words, \code{library(shinybootstrap2)}, shouldn't appear in your code.
#' This is because attaching the package will result in functions from
#' shinybootstrap2 masking functions from shiny, even outside of
#' \code{withBootstrap2()}.
#'
#' @param x An expression to evaluate with Bootstrap 2 components.
#' @param env The environment in which to evaluate \code{x}.
#' @param quoted Treat \code{x} as a quoted expression. If \code{FALSE} (the
#'   default) \code{x} will be treated as an unquoted expression. If
#'   \code{TRUE}, the code should be the output of a \code{\link{quote}()}.
#' @examples
#' \dontrun{
#' library(shiny)
#'
#' ## Single-file app using Bootstrap 2 ===========================
#' shinybootstrap2::withBootstrap2({
#'   shinyApp(
#'     ui = fluidPage(
#'       numericInput("n", "n", 1),
#'       plotOutput("plot")
#'     ),
#'     server = function(input, output) {
#'       output$plot <- renderPlot( plot(head(cars, input$n)) )
#'     }
#'   )
#' })
#'
#'
#' ## App with server.R and UI. R =================================
#' ## ui.R
#' shinybootstrap2::withBootstrap2({
#'   fluidPage(
#'     selectInput("ui", "Input type", choices = c("numeric", "slider")),
#'     uiOutput("n_ui"),
#'     plotOutput("plot")
#'   )
#' })
#'
#' ## server.R
#' # In server.R, it's only necessary to wrap code in withBoostrap2()
#' # when renderUI() is used.
#' shinybootstrap2::withBootstrap2({
#'   function(input, output) {
#'     output$n_ui <- renderUI({
#'       if (input$ui == "numeric")
#'         numericInput("n", "n", 1)
#'       else if (input$ui == "slider")
#'         sliderInput("n", "n", 1, 10, value = 1)
#'     })
#'     output$plot <- renderPlot( plot(head(cars, input$n)) )
#'   }
#' })
#'
#' }
#' @export
withBootstrap2 <- function(x, env = parent.frame(), quoted = FALSE) {
  if (!quoted) x <- substitute(x)

  # Copy everything from the shinybootstrap2 environment to a new environment
  # which is a child of the calling environment.
  bs2env <- list2env(bs2exports(), parent = env)

  eval(x, bs2env)
}


# Return a list of all exported objects from this package.
bs2exports <- function() {
  if (is.null(cache$exports)) {
    cache$exports <- mget(getNamespaceExports("shinybootstrap2"), asNamespace("shinybootstrap2"))
  }

  cache$exports
}

# A cache for the list of exported objects from this package.
cache <- new.env()


.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "You probably do not want to attach this package (with library() or require()).",
    " Instead, you should use shinybootstrap2::withBootstrap2().",
    " You can hide this message with suppressPackageStartupMessages()."
  )
}
