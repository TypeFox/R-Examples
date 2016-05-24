#############################################################################
# out_plain.R
#############################################################################

out.plain.init <- function()
{
    return(list(math = out.default.math,
                names = out.default.names,
                specialchar = out.default.specialchar,
                value = out.default.value,
                concat = out.default.concat,
                subscript = out.default.subscript,
                superscript = out.default.superscript,
                bracket = out.default.bracket,
                above = out.default.above,
                below = out.default.below))
}
