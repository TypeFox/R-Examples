"modelSetOR" <-
function(factoryName, overRelaxation)
#   Set over-relaxed updating
{
    name <- sQuote(factoryName)
    command <- paste("UpdaterMethods.SetFactory(", name, 
                     ") ;UpdaterMethods.OverRelaxationGuard;",
                     "UpdaterMethods.SetOverRelaxation(", 
                     overRelaxation,
                     ")", sep = "")
    .CmdInterpreter(command)
}
