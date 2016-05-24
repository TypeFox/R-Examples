#############################################################################
# out_markdown.R
#############################################################################

#' @include out_default.R
#' @include out_html.R
#' @include utils.R
NULL

out.markdown.init <- function()
{
    return(list(math = out.markdown.math,
                names = out.default.names,
                specialchar = out.markdown.specialchar,
                value = out.default.value,
                concat = out.default.concat,
                subscript = out.html.subscript,
                superscript = out.html.superscript,
                bracket = out.default.bracket,
                above = out.default.above,
                below = out.default.below))
}

out.markdown.math <- function(..., mmode)
{
    if (mmode)
        return(paste0("$", ..., "$"))
    else
        return(paste0(...))
}

out.markdown.specialchar <- function(x)
{
    myspec <- c("CHI" = "&chi;")
    
    return(utils.symbols.replace(x, replacements = myspec))
}
