AUEformat <-
function (ASC, UNI, EXP, format = "expression")
{
    if (format == "expression")
        return(EXP)
    if (format == "charstring")
        return(UNI)
    if (format == "ASCII")
        return(ASC)
}
