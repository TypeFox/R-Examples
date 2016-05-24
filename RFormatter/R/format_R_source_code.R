#' Formats R source code.
#'
#' \code{format_R_source_code} is very much based on
#' \href{http://yihui.name/formatR}{formatR}, but tries to improve it by
#' heuristics. For example, spaces can be forced around the division operator
#' \code{/}.
#'
#' @param source String with source code to format. This not the filename of the
#'     source file.
#' @param formatR_arguments List of arguments passed to
#'     \href{http://yihui.name/formatR}{formatR} via its
#'     \code{\link[formatR]{tidy_source}}.
#' @param remove_trailing_whitespace Boolean: should horizontal whitespace at
#'     the end of lines be removed?
#' @param spaced_operators Vector of strings with operators around which spaces
#'     are forced.
#'
#' @return String with formatted source code.
#'
#' @examples
#'
#' format_R_source_code("if (b) { f() }")
#' format_R_source_code("if (b) { f()\n\nf() }")
#' format_R_source_code("p = 2", list(arrow = TRUE, width.cutoff = 80))
#' format_R_source_code("e^x", spaced_operators = c("/"))
#'
#' \dontrun{
#' format_R_source_code("f()", text = NULL)
#' format_R_source_code("f()", output = TRUE)
#' }
#'
#' @seealso \code{\link[formatR]{tidy_source}}
#'
#' @export
format_R_source_code <- function(source, formatR_arguments = list(), remove_trailing_whitespace = TRUE,
    spaced_operators = c("/", "%%", "%/%", ":", "^")) {
    # Keep the position of the unnamed arguments.
    full_formatR_arguments <- c(formatR_arguments, text = source, output = FALSE)
    formatR_result <- do.call(formatR::tidy_source, full_formatR_arguments)

    source_lines <- unlist(lapply(formatR_result$text.tidy, split_separated_lines))
    source_lines <- improve_formatR(source_lines, remove_trailing_whitespace, spaced_operators)

    newline <- guess_default_newline(source)
    paste(source_lines, collapse = newline)
}

improve_formatR <- function(source_lines, remove_trailing_whitespace, spaced_operators) {
    # Remove horizontal whitespace at the end of lines.
    if (remove_trailing_whitespace) {
        trailing_whitespace <- "\\h+$"
        source_lines <- replace_pattern(trailing_whitespace, "", source_lines)
    }

    source_lines <- space_operators(source_lines, spaced_operators)

    source_lines
}


# Tries to force spaces around certain binary operators.
space_operators <- function(source_lines, spaced_operators) {
    if (length(spaced_operators) != 0) {
        # Regular expression for leaves of the abstract syntax tree.
        leaf_character <- "[\\w.]"
        before <- paste0(leaf_character, "|[])]")
        spaced_operator <- paste(sapply(spaced_operators, make_regular_expression_literal),
            collapse = "|")
        after <- paste0(leaf_character, "|[(+-]")
        search <- paste0("(", before, ")(", spaced_operator, ")(", after, ")")
        replacement <- "\\1 \\2 \\3"

        space_actual_code <- function(actual_code) {
            replace_pattern(search, replacement, actual_code)
        }

        source_lines <- sapply(source_lines, function(source_line) {
            transform_R_source_line(space_actual_code, source_line)
        })
    }

    source_lines
}
