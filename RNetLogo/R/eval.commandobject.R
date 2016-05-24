eval.commandobject <- function(command) {
    # if the command is of length 1, just return it
    # else transform it to a NetLogo-List
    if ((is.vector(command)) & (length(command) > 1))
    {      
      # recursive call (needed for lists), create nested NetLogo-Lists
      if (is.list(command))
      {
        command <- lapply(command, function(x) {eval.commandobject(x)})
      }
      # if the elements are of type character, then create the quotes around each element
      if (is.character(command))
      {
        command <- paste(command, collapse="\" \"")
        command <- paste(c('[ ','\"',command,'\"',' ]'), collapse="")
      }
      # if the elements aren't character just paste them together
      else
      {
        command <- paste(c('[',command,']'), collapse=" ")
      }
    }
    # if the command is a dataframe
    if (is.data.frame(command))
    {
      # recursive call for each column of the dataframe (creates nested NetLogo-Lists)
      # TODO: how to handel level variables??
      command <- lapply(command, function(x) {eval.commandobject(x)})
      # at the end, paste base brackets around the nested NetLogo-List
      command <- paste(c('[',command,']'), collapse=" ")
    }
    # if the command is a level/factor (in dataframe) convert it into character 
    if (is.factor(command))
    {
      command <- eval.commandobject(as.character(command))
    }
    return (command)
}