"modelNames" <-
function()
{
#   gets names in OpenBUGS model
    command <- "BugsRobjects.GetNumberNames"
    number <- .Integer(command)
    name <- character(number)
    if(length(number)){
      cmds <- character(0)
      cmdtype <- character()
      for(i in 1:number){
        cmds <- c(cmds, paste("BugsRobjects.SetIndex(", i-1, ")", sep=""),
                  "BugsRobjects.GetStringLength")
        cmdtype <- c(cmdtype, c("CmdInterpreter","Integer"))
      }    
      res <- .OpenBUGS(cmds, cmdtype)
      numchar <- unlist(res[seq(2, 2*number, by=2)])
        
      cmds <- character(0)
      cmdtype <- character()
      args <- list()
      for(i in 1:number){
        char <- paste(rep(" ", numchar[i]), collapse="")
        cmds <- c(cmds,
                  paste("BugsRobjects.SetIndex(", i-1, ")", sep=""),
                  "BugsRobjects.GetVariable")
        cmdtype <- c(cmdtype, c("CmdInterpreter","CharArray"))
        args <- c(args, list(NA, char))
      }       
      res <- .OpenBUGS(cmds, cmdtype, args)
      name <- unlist(res[seq(2, 2*number, by=2)])      
    }
    return(name)
  }
