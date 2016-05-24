"modelPrecision" <-
function(prec)
#   Set the precision to which results are displayed
{
    if(!is.numeric(prec))
        stop("prec ", "must be numeric")
    prec <- as.integer(prec)
    options(BRugsPrec=prec)
#    command <- paste("BugsMappers.SetPrec(", prec, ")")
#    invisible(.CmdInterpreter(command))
}
