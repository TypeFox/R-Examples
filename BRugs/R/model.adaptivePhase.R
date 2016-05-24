"modelAdaptivePhase" <-
function()
#   Get endOfAdapting field
{
    command <- "BugsInterface.endOfAdapting"
    .Integer(command) - 1
}
