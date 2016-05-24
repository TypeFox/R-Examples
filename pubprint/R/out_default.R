#############################################################################
# out_default.R
#############################################################################

#' @include utils.R
NULL

out.default.names <- function(x)
{
    ret <- c()

    for (xi in x)
        ret <- c(ret,
                 switch(EXPR = xi,
                        "cor" = "r", 
                        "D^-" = paste0("D", out.superscript("-")),
                        "X-squared" = paste0(out.specialchar("CHI"),
                                             out.superscript("2")), 
                        "Bartlett's K-squared" = paste0(out.specialchar("CHI"),
                                                        out.superscript("2")), 
                        "Friedman chi-squared" = paste0(out.specialchar("CHI"),
                                                        out.superscript("2")),
                        "mean of x" = paste0("M", out.subscript("x")), 
                        "mean of y" = paste0("M", out.subscript("y")),
                        xi))

    return(ret)
}

out.default.specialchar <- function(x)
{
    return(x)
}

out.default.math <- function(..., mmode)
{
    if (missing(mmode))
        stop("argument \"mmode\" is missing, with no default")

    paste0(...)
}

out.default.value <- function(x,
                              name,
                              inbracket,
                              nsmall,
                              replace0,
                              leading0,
                              drop0trailing)
{
    if (!missing(name) && length(x) > 1 && length(name) > 1 && length(x) !=
        length(name))
        stop("argument name must be length one or length of argument x")
    if (!missing(inbracket) && missing(name))
        stop("If argument inbracket is given, argument name must be given as well")
    if (!missing(inbracket) && length(inbracket) != length(name) && length(inbracket) > 1)
        stop("argument inbracket must be length one or equal to length of argument name")

    num <- x

    if (replace0)
    {
        num <- ifelse(abs(x) < 1/(10^nsmall), 
                      1/(10^nsmall) * sign(x),
                      num)
        # cases with x = 0 we have to treat specially
        num <- ifelse(num != 0, 
                      num,
                      1/(10^nsmall))
    }

    num <- format(round(num, nsmall), 
                  nsmall = nsmall, 
                  trim = TRUE,
                  drop0trailing = drop0trailing)

    if (!leading0) 
        num <- sub("^(-?)0\\.", "\\1\\.", num)

    if (!missing(inbracket))
        name <- paste0(name, out.bracket(inbracket))

    if (missing(name) || is.null(name))
        num <- paste(num)
    else
        num <- ifelse(abs(x) < 1/(10^nsmall),
                      paste(name, num, sep = out.specialchar("<")),
                      paste(name, num, sep = out.specialchar("=")))

    return(num)
}

out.default.concat <- function(..., sep)
{
    stringr::str_c(..., sep = sep, collapse = sep)
}

out.default.subscript <- function(x)
{
    return(x)
}

out.default.superscript <- function(x)
{
    return(x)
}

out.default.bracket <- function(x, brackets)
{
    sapply(x, 
           out.default.bracket.work, brackets = brackets, 
           USE.NAMES = FALSE)
}

out.default.bracket.work <- function(x, brackets)
{
    if (length(brackets) %% 2 != 0 && length(brackets) != 1)
        stop("Argument brackets must be length one or a multiple of two.")
    if (1 != length(x))
        stop("Argument x must be length one.")

    if (1 == length(brackets))
        return(paste0(brackets, x, brackets))

    positions <- stringr::str_locate_all(x, stringr::fixed(brackets))

    for (i in seq_len(length(brackets)))
    {
        bracket_new <- brackets[(i + 1) %% length(brackets) + 1]

        mypos <- positions[[i]][, "start"]

        for (j in seq_len(length(mypos)))
            stringr::str_sub(x, start = mypos[j], end = mypos[j]) <- bracket_new
    }

    return(paste0(brackets[1], x, brackets[2]))
}

out.default.above <- function(x, y)
{
    paste0(x, y)
}

out.default.below <- function(x, y)
{
    paste0(x, y)
}
