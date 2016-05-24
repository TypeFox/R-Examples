"modelSetIts" <-
function(factoryName, iterations)
#   Set maximum number of iterations in iterative algorithms
{
    name <- sQuote(factoryName)
    command <- paste("UpdaterMethods.SetFactory(", name, 
                     ") ;UpdaterMethods.IterationsGuard;",
                     "UpdaterMethods.SetIterations(", 
                     iterations,
                     ")", sep = "")
    .CmdInterpreter(command)
}
