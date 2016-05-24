# Shiny objects if the widget sets elementId, so we
# need to detect it.  Thanks to Joe Cheng for suggesting this code.
inShiny <- function() !is.null(getDefaultReactiveDomain())

rglwidget <- local({
  reuseDF <- NULL

  function(x = scene3d(), width = NULL, height = NULL,
           controllers = NULL, snapshot = FALSE,
           elementId = NULL,
           reuse = !interactive(), ...) {
  if (is.na(reuse))
    reuseDF <- NULL # local change only
  else if (!reuse)
    reuseDF <<- NULL

  if (is.null(elementId) && !inShiny())
    elementId <- paste0("rgl", sample(100000, 1))

  x = convertScene(x, width, height, snapshot = snapshot,
                   elementId = elementId, reuse = reuseDF)
  if (!is.na(reuse))
    reuseDF <<- attr(x, "reuse")

  if (inherits(controllers, "rglPlayer") && is.na(controllers$x$sceneId)) {
    x$players <- controllers$elementId
    if (!is.null(elementId))
      controllers$x$sceneId <- elementId
  } else if (inherits(controllers, "shiny.tag.list")) {
    x$players <- character()
    for (i in seq_along(controllers)) {
      if (inherits(controllers[[i]], "rglPlayer") && is.na(controllers[[i]]$x$sceneId)) {
        x$players <- c(x$players, controllers[[i]]$elementId)
        if (!is.null(elementId))
          controllers[[i]]$x$sceneId <- elementId
      }
    }
  } else if (!is.null(controllers))
    x$players <- controllers
  # create widget
  result <- structure(htmlwidgets::createWidget(
    name = 'rglWebGL',
    x = x,
    width = width,
    height = height,
    package = 'rglwidget',
    elementId = elementId,
    ...
  ), rglReuse = attr(x, "reuse"))
  if (inherits(controllers, "shiny.tag.list")) {
    controllers[[length(controllers) + 1]] <- result
    controllers
  } else if (inherits(controllers, "rglPlayer")) {
    browsable(tagList(controllers, result))
  } else
    result
  }
})

#' Widget output function for use in Shiny
#'
#' @export
rglwidgetOutput <- function(outputId, width = '512px', height = '512px'){
  shinyWidgetOutput(outputId, 'rglWebGL', width, height, package = 'rglwidget')
}

#' Widget render function for use in Shiny
#'
#' @export
renderRglwidget <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, rglwidgetOutput, env, quoted = TRUE)
}
