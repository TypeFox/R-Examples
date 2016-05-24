#############################################################################
# out_latex.R
#############################################################################

#' @include out_default.R
#' @include utils.R
NULL

out.latex.init <- function()
{
    return(list(math = out.latex.math,
                names = out.default.names,
                specialchar = out.latex.specialchar,
                value = out.default.value,
                concat = out.default.concat,
                subscript = out.latex.subscript,
                superscript = out.latex.superscript,
                bracket = out.default.bracket,
                above = out.latex.above,
                below = out.latex.below))
}

out.latex.specialchar <- function(x)
{
    myspec <- c("<" = "\\ifmmode<\\else\\textless\\fi",
                ">" = "\\ifmmode<\\else\\textgreater\\fi",
                "CHI" = "\\ifmmode\\chi\\else\\(\\chi\\)\\fi")

    return(utils.symbols.replace(x, replacements = myspec))
}

out.latex.math <- function(..., mmode)
{
    if (mmode)
        return(paste0("\\ensuremath{", ..., "}"))
    else
        return(paste0(...))
}

out.latex.subscript <- function(x)
{
    return(paste0("\\ifmmode_{", 
                  x, 
                  "}\\else\\textsubscript{", 
                  x,
                  "}\\fi"))
}

out.latex.superscript <- function(x)
{
    return(paste0("\\ifmmode^{", 
                  x, 
                  "}\\else\\textsuperscript{", 
                  x,
                  "}\\fi"))
}

# requires amsmath package
# x over y
out.latex.above <- function(x, y)
{
    if ("^" == y)
        return(paste0("\\ifmmode\\hat{",
                      x, 
                      "}\\else ",
                      x,
                      "\\textsuperscript{\\textasciicircum}\\fi"))
    else 
        return(paste0("\\ifmmode\\overset{",
                      y, 
                      "}{", 
                      x,
                      "}\\else ",
                      x,
                      "\\textsuperscript{",
                      y,
                      "}\\fi"))
}

# requires amsmath package
# x below y
out.latex.below <- function(x, y)
{
    return(paste0("\\ifmmode\\underset{",
                  y, 
                  "}{", 
                  x,
                  "}\\else ",
                  x,
                  "\\textsubscript{",
                  y,
                  "}\\fi"))
}
