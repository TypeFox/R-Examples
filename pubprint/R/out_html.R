#############################################################################
# out_html.R
#############################################################################

#' @include out_default.R
NULL

out.html.init <- function()
{
    return(list(math = out.html.math,
                names = out.default.names,
                specialchar = out.html.specialchar,
                value = out.default.value,
                concat = out.default.concat,
                subscript = out.html.subscript,
                superscript = out.html.superscript,
                above = out.default.above,
                below = out.default.below))
}

out.html.specialchar <- function(x)
{
    myspec <- c("<" = "&lt;",
                ">" = "&gt;",
                "CHI" = "&chi;")
    
    return(utils.symbols.replace(x, replacements = myspec))
}

out.html.math <- function(..., mmode)
{
    if (mmode)
        return(paste0("<MATH>", ..., "</MATH>"))
    else
        return(paste0(...))
}

out.html.subscript <- function(x)
{
    return(paste0("<sub>", x, "</sub>"))
}

out.html.superscript <- function(x)
{
    return(paste0("<sup>", x, "</sup>"))
}
