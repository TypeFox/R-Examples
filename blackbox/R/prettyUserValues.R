prettyUserValues <-
function (innervalues, extradigits = 0)
{
    uservalues <- toUserValues(innervalues)
    return(prettynum(uservalues, extradigits = extradigits))
}
