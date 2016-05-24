"modelSetAP" <-
function(factoryName, adaptivePhase)
#   Set the length of adaptive phase
{
    name <- sQuote(factoryName)
    command <- paste("UpdaterMethods.SetFactory(", name, 
                     ") ;UpdaterMethods.AdaptivePhaseGuard;",
                     "UpdaterMethods.SetAdaptivePhase(", 
                     adaptivePhase,
                     ")", sep = "")
    .CmdInterpreter(command)
}
