Tree <-
function(value, children = list())
{
    y <- list(value = value, children = as.list(children))
    class(y) <- "Tree"
    y
}

format.Tree <-
function(x,
         width = 0.9 * getOption("width"), indent = 0,
         brackets = c("(", ")"), ...)
{
    ffmt <- function(x) {
        sprintf("%s%s %s%s",
                brackets[1L],
                x$value,
                paste(sapply(x$children,
                             function(e) {
                                 if(inherits(e, "Tree"))
                                     ffmt(e)
                                 else
                                     format(e)
                             }),
                      collapse = " "),
                brackets[2L])
    }

    s <- ffmt(x)
    if(nchar(s) + indent < width) return(s)

    y <- sapply(x$children,
                function(e) {
                    if(inherits(e, "Tree"))
                        format(e, width = width, indent = indent + 2L,
                               brackets = brackets)
                    else
                        format(e)
                })
    y <- sprintf("\n%s%s",
                 paste(rep.int(" ", indent + 2L), collapse = ""),
                 y)
    sprintf("%s%s%s%s",
            brackets[1L],
            x$value,
            paste(y, collapse = ""),
            brackets[2L])
}

## print.Tree <-
## function(x, ...)
## {
##     writeLines(format(x, ...))
##     invisible(x)
## }
    
Tree_parse <-
function(x, brackets = c("(", ")"))
{
    errfmt <- function(token, expected) {
        sprintf("expected %s but got %s", expected, token)
    }

    re_o <- sprintf("\\%s", brackets[1L]) # open
    re_c <- sprintf("\\%s", brackets[2L]) # close
    re_n <- sprintf("[^\\s%s%s]+", re_o, re_c) # node
    re_l <- sprintf("[^\\s%s%s]+", re_o, re_c) # leaf
    re <- sprintf("%s\\s*(%s)?|%s|(%s)", re_o, re_n, re_c, re_l)
    m <- gregexpr(re, x, perl = TRUE)

    stack <- list(list(NULL, list()))

    for(token in regmatches(x, m)[[1L]]) {
        if(substring(token, 1L, 1L) == "(") {
            if((length(stack) == 1L) &&
               (length(stack[[1L]][[2L]]) > 0L))
                stop(errfmt(sQuote(token), "end of string"))
            value <- sub("\\s*", "", substring(token, 2L))
            stack <- c(stack, list(Tree(value, list())))
        }
        else if(token == ")") {
            if((n <- length(stack)) == 1L) {
                if(!length(stack[[1L]][[2L]]))
                    stop(errfmt(sQuote(token), sQuote(brackets[1L])))
                else
                    stop(errfmt(sQuote(token), "end of string"))
            }
            elt <- stack[[n]]
            ## class(elt) <- "Tree"
            stack <- stack[-n]
            n <- n - 1L
            stack[[n]][[2L]] <- c(stack[[n]][[2L]], list(elt))
        }
        else {
            if((n <- length(stack)) == 1L)
                stop(errfmt(sQuote(token), sQuote(brackets[1L])))
            stack[[n]][[2L]] <- c(stack[[n]][[2L]], list(token))
        }
    }

    if(length(stack) > 1L)
        stop(errfmt("end of string", sQuote(brackets[2L])))
    else if(!length(stack[[1L]][[2L]]))
        stop(errfmt("end of string", sQuote(brackets[1L])))

    stack[[1L]][[2L]][[1L]]
}

Tree_apply <-
function(x, f, recursive = FALSE)
{
    if(!recursive) return(lapply(x$children, f))

    g <- function(e) {
        y <- f(e)
        if(inherits(e, "Tree"))
            list(y, lapply(e$children, g))
        else
            y
    }

    lapply(x$children, g)
}
