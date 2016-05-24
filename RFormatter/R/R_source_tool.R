# Replaces actual code in a line of R source code.
transform_R_source_line <- function(transform_actual_code, source_line) {
    change <- 0
    change_indices <- find_actual_R_code(source_line)
    change_count <- length(change_indices)
    transformed_line <- ""

    while (change <= change_count) {
        if (change >= 1) {
            start <- change_indices[change]
        } else {
            start <- 1
        }
        if (change < change_count) {
            end <- change_indices[change + 1] - 1
        } else {
            end <- nchar(source_line)
        }

        source_part <- substr(source_line, start, end)
        if (change %% 2 == 0) {
            # Actual code to replace.
            source_part <- transform_actual_code(source_part)
        }
        transformed_line <- paste0(transformed_line, source_part)

        change <- change + 1
    }

    transformed_line
}

# Tries to distinguish quotes and comments from any other kind of code.
find_actual_R_code <- function(source_line) {
    quotes <- c("\"", "'", "`")
    signed_changes <- c()

    # Register index at a border.
    add_change <- function(index, is_end) {
        if (is_end) {
            # End.
            signed_index <- -index
        } else {
            # Start.
            signed_index <- index
        }
        c(signed_changes, signed_index)
    }

    # Iteration over the characters.
    index <- 1
    count <- nchar(source_line)

    # In which quote are we, if any at all?
    maybe_quote_index <- NA
    # If we are in a quote, is the current character escaped?
    is_escaped <- FALSE
    # Are we in a comment?
    is_comment <- FALSE

    # We can terminate early if we see a line comment.
    while (index <= count && !is_comment) {
        character <- substr(source_line, index, index)

        if (is.na(maybe_quote_index)) {
            # Not in a quote.
            maybe_quote_index <- match(character, quotes)
            if (!is.na(maybe_quote_index)) {
                # Start of a quote.
                signed_changes <- add_change(index, FALSE)
            }
            if (character == "#") {
                # Start of a comment.
                signed_changes <- add_change(index, FALSE)
                is_comment <- TRUE
            }
        } else if (is_escaped) {
            # Escaped character in a quote.
            is_escaped <- FALSE
        } else {
            # In a quote without escaping.
            quote_now <- match(character, quotes)
            if (!is.na(quote_now) && quote_now == maybe_quote_index) {
                # End of a quote.
                signed_changes <- add_change(index, TRUE)
                maybe_quote_index <- NA
            }
            is_escaped <- character == "\\"
        }

        index <- index + 1
    }

    make_change_indices_unsigned(signed_changes)
}

# Simplifies start and end indices to start indices only.
make_change_indices_unsigned <- function(signed_changes) {
    unsigned_changes <- c()

    index <- 1
    count <- length(signed_changes)

    while (index <= count) {
        signed_change <- signed_changes[index]

        if (signed_change > 0) {
            # Start.
            changes <- signed_change
            index_gap <- 1
        } else {
            # End.
            unsigned_change <- -signed_change + 1
            if (index < count && signed_changes[index + 1] == unsigned_change) {
                # End is immediately followed by a start: remove both.
                changes <- c()
                index_gap <- 2
            } else {
                # Replace the end here by a start at next index.
                changes <- unsigned_change
                index_gap <- 1
            }
        }

        unsigned_changes <- c(unsigned_changes, changes)
        index <- index + index_gap
    }

    unsigned_changes
}
