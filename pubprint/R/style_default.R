#############################################################################
# style_default.R
#############################################################################

#########
# initialise function
#########

# not necessary

#########
# working functions
#########

# character function
# x: assuming all list items are characters
# return: list items as a vector
style.default.character <- function(x)
{
    return(unlist(x))
}

# numeric function
# x: assuming all list items are numeric
# return: vector from out.value
style.default.numeric <- function(x, ...)
{
    return(out.value(unlist(x), ...))
}
