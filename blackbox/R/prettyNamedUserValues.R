prettyNamedUserValues <-
function (innervalues, extradigits = 0)
{
    uservalues <- prettyUserValues(innervalues, extradigits = extradigits)
    names(uservalues) <- sapply(names(innervalues), formatName,
        format = "ASCII")
    return(uservalues)
}
