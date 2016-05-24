modelDisable <- function(factory){
    command <- paste("UpdaterMethods.SetFactory('", factory,"');UpdaterMethods.Disable", sep = "")
    invisible(.CmdInterpreter(command))

}


modelEnable <- function(factory){
    command <- paste("UpdaterMethods.SetFactory('", factory,"');UpdaterMethods.Enable", sep = "")
    invisible(.CmdInterpreter(command))
}
