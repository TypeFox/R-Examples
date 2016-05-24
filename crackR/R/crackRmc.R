crackRmc <-
function(parameters)
{
    analysis.type <- parameters$analysis.type

    if( analysis.type == "single" )
        {
            class(parameters) <- "Sing"
        }

    if( analysis.type == "multiple" )
        {
            class(parameters) <- "Mult"
        }

    if( analysis.type == "CD" )
        {
            class(parameters) <- "CD"
        }

    UseMethod("crackRmc", parameters)
}
