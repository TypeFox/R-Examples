most_frequent_value <- function(values) {
    names(which.max(table(values)))
}

make_regular_expression_literal <- function(string) {
    character_to_escape <- "([$()*+.?[\\^{|])"
    escaped_character <- "\\\\\\1"
    regular_expression <- replace_pattern(character_to_escape, escaped_character,
        string)
}

replace_pattern <- function(search, replacement, string) {
    gsub(search, replacement, string, perl = TRUE)
}

# Breaks a string up into its lines. Any of LF, CR, and CR+LF are treated as
# newlines. Unlike \code{\link[base]{readLines}}, this interprets a newline as
# a separator, not a terminator.
split_separated_lines <- function(string) {
    scan(what = character(), sep = "\n", na.strings = c(), blank.lines.skip = FALSE,
        quiet = TRUE, text = string)
}

newline_pattern <- function() {
    # Do not use the vertical whitespace "\\v" here, since it treats "\r\n" as
    # two instead of a single newline. For the same reason, the order of the
    # alternatives in the following regular expression matters.
    "\\r?\\n|\\r"
}

guess_default_newline <- function(string) {
    matches <- gregexpr(newline_pattern(), string, perl = TRUE)
    if (matches[[1]][1] == -1) {
        # Default.
        "\n"
    } else {
        newlines <- regmatches(string, matches)
        most_frequent_value(unlist(newlines))
    }
}
