#' Create a select list input control
#'
#' Create a select list that can be used to choose a single or
#' multiple items from a list of values.
#'
#' By default, \code{selectInput()} and \code{selectizeInput()} use the
#' JavaScript library \pkg{selectize.js} (\url{https://github.com/brianreavis/selectize.js})
#' to instead of the basic select input element. To use the standard HTML select
#' input element, use \code{selectInput()} with \code{selectize=FALSE}.
#'
#' @param inputId Input variable to assign the control's value to
#' @param label Display label for the control, or \code{NULL}
#' @param choices List of values to select from. If elements of the list are
#' named then that name rather than the value is displayed to the user.
#' @param selected The initially selected value (or multiple values if
#' \code{multiple = TRUE}). If not specified then defaults to the first value
#' for single-select lists and no values for multiple select lists.
#' @param multiple Is selection of multiple items allowed?
#' @param selectize Whether to use \pkg{selectize.js} or not.
#' @return A select list control that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link[shiny]{updateSelectInput}}
#'
#' @examples
#' selectInput("variable", "Variable:",
#'             c("Cylinders" = "cyl",
#'               "Transmission" = "am",
#'               "Gears" = "gear"))
#' @export
selectInput <- function(inputId, label, choices, selected = NULL,
                        multiple = FALSE, selectize = TRUE, width = NULL) {
  # resolve names
  choices <- choicesWithNames(choices)

  # default value if it's not specified
  if (is.null(selected)) {
    if (!multiple) selected <- firstChoice(choices)
  } else selected <- validateSelected(selected, choices, inputId)

  # create select tag and add options
  selectTag <- tags$select(id = inputId, selectOptions(choices, selected))
  if (multiple)
    selectTag$attribs$multiple <- "multiple"

  # return label and select tag
  res <- tagList(controlLabel(inputId, label), selectTag)
  if (!selectize) return(res)
  selectizeIt(inputId, res, NULL, width, nonempty = !multiple && !("" %in% choices))
}

firstChoice <- function(choices) {
  if (length(choices) == 0L) return()
  choice <- choices[[1]]
  if (is.list(choice)) firstChoice(choice) else choice
}

# Create tags for each of the options; use <optgroup> if necessary.
# This returns a HTML string instead of tags, because of the 'selected'
# attribute.
selectOptions <- function(choices, selected = NULL) {
  html <- mapply(choices, names(choices), FUN = function(choice, label) {
    if (is.list(choice)) {
      # If sub-list, create an optgroup and recurse into the sublist
      sprintf(
        '<optgroup label="%s">\n%s\n</optgroup>',
        htmlEscape(label),
        selectOptions(choice, selected)
      )

    } else {
      # If single item, just return option string
      sprintf(
        '<option value="%s"%s>%s</option>',
        htmlEscape(choice),
        if (choice %in% selected) ' selected' else '',
        htmlEscape(label)
      )
    }
  })

  HTML(paste(html, collapse = '\n'))
}

#' @rdname selectInput
#' @param ... Arguments passed to \code{selectInput()}.
#' @param options A list of options. See the documentation of \pkg{selectize.js}
#'   for possible options (character option values inside \code{\link{I}()} will
#'   be treated as literal JavaScript code; see \code{\link[shiny]{renderDataTable}()}
#'   for details).
#' @param width The width of the input, e.g. \code{'400px'}, or \code{'100\%'};
#'   see \code{\link{validateCssUnit}}.
#' @note The selectize input created from \code{selectizeInput()} allows
#'   deletion of the selected option even in a single select input, which will
#'   return an empty string as its value. This is the default behavior of
#'   \pkg{selectize.js}. However, the selectize input created from
#'   \code{selectInput(..., selectize = TRUE)} will ignore the empty string
#'   value when it is a single choice input and the empty string is not in the
#'   \code{choices} argument. This is to keep compatibility with
#'   \code{selectInput(..., selectize = FALSE)}.
#' @export
selectizeInput <- function(inputId, ..., options = NULL, width = NULL) {
  selectizeIt(inputId, selectInput(inputId, ..., selectize = FALSE), options, width)
}

# given a select input and its id, selectize it
selectizeIt <- function(inputId, select, options, width = NULL, nonempty = FALSE) {
  res <- checkAsIs(options)

  selectizeDep <- htmlDependency(
    "selectize", "0.8.5",
    c(file = system.file("www/selectize", package = "shinybootstrap2")),
    stylesheet = "css/selectize.bootstrap2.css",
    head = format(tagList(
      HTML('<!--[if lt IE 9]>'),
      tags$script(src = 'selectize-0.8.5/js/es5-shim.min.js'),
      HTML('<![endif]-->'),
      tags$script(src = 'selectize-0.8.5/js/selectize.min.js')
    ))
  )
  attachDependencies(
    tagList(
      select,
      tags$script(
        type = 'application/json',
        `data-for` = inputId, `data-nonempty` = if (nonempty) '',
        `data-eval` = if (length(res$eval)) HTML(toJSON(res$eval)),
        `data-width` = validateCssUnit(width),
        if (length(res$options)) HTML(toJSON(res$options)) else '{}'
      )
    ),
    selectizeDep
  )
}

# for options passed to DataTables/Selectize/..., the options of the class AsIs
# will be evaluated as literal JavaScript code
checkAsIs <- function(options) {
  evalOptions <- if (length(options)) {
    nms <- names(options)
    if (length(nms) == 0L || any(nms == '')) stop("'options' must be a named list")
    i <- unlist(lapply(options, function(x) {
      is.character(x) && inherits(x, 'AsIs')
    }))
    if (any(i)) {
      # must convert to character, otherwise toJSON() turns it to an array []
      options[i] <- lapply(options[i], paste, collapse = '\n')
      nms[i]  # options of these names will be evaluated in JS
    }
  }
  list(options = options, eval = evalOptions)
}
