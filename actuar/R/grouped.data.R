### ===== actuar: An R Package for Actuarial Science =====
###
### Creation of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon, Louis-Philippe Pouliot

grouped.data <- function(..., right = TRUE, row.names = NULL,
                         check.rows = FALSE, check.names = TRUE)
{
    ## Utility function
    numform <- function(x, w)
        formatC(x, digits = 2, width = w, format = "fg")

    ## The function must be called with at least two arguments. The
    ## first is the vector of group boundaries. The others are vectors
    ## of group frequencies. All arguments will be converted to data
    ## frames.
    x <- list(...)
    xnames <- names(x)                  # preserve names
    y <- as.data.frame(x[-1])           # group frequencies
    x <- as.data.frame(x[[1]])          # group boundaries
    nx <- nrow(x)
    ny <- nrow(y)

    ## There must be exactly one group boundary more than frequencies.
    if (nx - ny != 1)
        stop("invalid number of group boundaries and frequencies")

    ## Replace missing frequencies by zeros.
    nax <- is.na(x)
    if (any(nax))
    {
        x[nax] <- 0
        warning("missing frequencies replaced by zeros")
    }

    ## Return a data frame with formatted group boundaries in the
    ## first column.
    w <- max(nchar(x[-1, ]))            # longest upper boundary
    xfmt <- paste(if (right) "(" else "[",
                  numform(x[-nx, ], -1), ", ", numform(x[-1, ], w),
                  if (right) "]" else ")",
                  sep = "")
    res <- data.frame(xfmt, y, row.names = row.names, check.rows = check.rows,
                      check.names = check.names)
    names(res) <- c(xnames[1], names(y))
    class(res) <- c("grouped.data", "data.frame")
    environment(res) <- new.env()
    assign("cj", unlist(x, use.names = FALSE), environment(res))
    attr(res, "right") <- right
    res
}
