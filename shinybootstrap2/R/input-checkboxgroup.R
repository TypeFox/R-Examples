#' Checkbox Group Input Control
#'
#' Create a group of checkboxes that can be used to toggle multiple choices
#' independently. The server will receive the input as a character vector of the
#' selected values.
#'
#' @param inputId Input variable to assign the control's value to.
#' @param label Display label for the control, or \code{NULL}.
#' @param choices List of values to show checkboxes for. If elements of the list
#'   are named then that name rather than the value is displayed to the user.
#' @param selected The values that should be initially selected, if any.
#' @param inline If \code{TRUE}, render the choices inline (i.e. horizontally)
#' @return A list of HTML elements that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link{checkboxInput}}, \code{\link{updateCheckboxGroupInput}}
#'
#' @examples
#' checkboxGroupInput("variable", "Variable:",
#'                    c("Cylinders" = "cyl",
#'                      "Transmission" = "am",
#'                      "Gears" = "gear"))
#'
#' @export
checkboxGroupInput <- function(inputId, label, choices, selected = NULL, inline = FALSE) {
  # resolve names
  choices <- choicesWithNames(choices)
  if (!is.null(selected))
    selected <- validateSelected(selected, choices, inputId)

  options <- generateOptions(inputId, choices, selected, inline)

  # return label and select tag
  tags$div(id = inputId,
           class = "control-group shiny-input-checkboxgroup",
           controlLabel(inputId, label),
           options)
}


# need <optgroup> when choices contains sub-lists
needOptgroup <- function(choices) {
  any(vapply(choices, is.list, logical(1)))
}

# Before shiny 0.9, `selected` refers to names/labels of `choices`; now it
# refers to values. Below is a function for backward compatibility.
validateSelected <- function(selected, choices, inputId) {
  # drop names, otherwise toJSON() keeps them too
  selected <- unname(selected)
  # if you are using optgroups, you're using shiny > 0.10.0, and you should
  # already know that `selected` must be a value instead of a label
  if (needOptgroup(choices)) return(selected)

  if (is.list(choices)) choices <- unlist(choices)

  nms <- names(choices)
  # labels and values are identical, no need to validate
  if (identical(nms, unname(choices))) return(selected)
  # when selected labels instead of values
  i <- (selected %in% nms) & !(selected %in% choices)
  if (any(i)) {
    warnFun <- if (all(i)) {
      # replace names with values
      selected <- unname(choices[selected])
      warning
    } else stop  # stop when it is ambiguous (some labels == values)
    warnFun("'selected' must be the values instead of names of 'choices' ",
            "for the input '", inputId, "'")
  }
  selected
}

# generate options for radio buttons and checkbox groups (type = 'checkbox' or
# 'radio')
generateOptions <- function(inputId, choices, selected, inline, type = 'checkbox') {
  # create tags for each of the options
  ids <- paste0(inputId, seq_along(choices))
  # generate a list of <input type=? [checked] />
  mapply(
    ids, choices, names(choices),
    FUN = function(id, value, name) {
      inputTag <- tags$input(
        type = type, name = inputId, id = id, value = value
      )
      if (value %in% selected)
        inputTag$attribs$checked <- "checked"
      tags$label(
        class = paste(type, if (inline) "inline"),
        inputTag, tags$span(name)
      )
    },
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
}

# Takes a vector or list, and adds names (same as the value) to any entries
# without names.
choicesWithNames <- function(choices) {
  # Take a vector or list, and convert to list. Also, if any children are
  # vectors with length > 1, convert those to list. If the list is unnamed,
  # convert it to a named list with blank names.
  listify <- function(obj) {
    # If a list/vector is unnamed, give it blank names
    makeNamed <- function(x) {
      if (is.null(names(x))) names(x) <- character(length(x))
      x
    }

    res <- lapply(obj, function(val) {
      if (is.list(val))
        listify(val)
      else if (length(val) == 1 && is.null(names(val)))
        val
      else
        makeNamed(as.list(val))
    })

    makeNamed(res)
  }

  choices <- listify(choices)
  if (length(choices) == 0) return(choices)

  # Recurse into any subgroups
  choices <- mapply(choices, names(choices), FUN = function(choice, name) {
    if (!is.list(choice)) return(choice)
    if (name == "") stop('All sub-lists in "choices" must be named.')
    choicesWithNames(choice)
  }, SIMPLIFY = FALSE)

  # default missing names to choice values
  missing <- names(choices) == ""
  names(choices)[missing] <- as.character(choices)[missing]

  choices
}
