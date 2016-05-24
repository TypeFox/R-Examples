NLSetAgentSet <-
function(agentset, input, var.name=NULL, nl.obj=NULL) 
{
  # get internal nl.obj if NULL
  if (is.null(nl.obj))
  {
    nl.obj <- "_nl.intern_"
  }
  # check for unknown nl.obj
  if (!(nl.obj %in% .rnetlogo$objects)) {
    stop(paste('There is no NetLogo reference stored under the name ',nl.obj,".", sep=""))
  }  
  
  if (is.data.frame(input)) {    
    agentset.len = NLReport(paste("count ",agentset,sep=" "), nl.obj=nl.obj)    
    if (length(input[[1]]) != agentset.len) {
      stop("Length of agentset not equal to length of input.")
    }    
    # get agent variable names
    vars_ <- names(input)
    prev_ <- paste("(foreach sort ",agentset, " ", sep="")
    ask_ <- "[ask ?1 ["
    cnt = 1
    inp_ <- paste(lapply(input,
                         function(x) 
                           paste("[",paste(x, collapse=" "), "]", collapse=" ")
                  ), " ", collapse=" ")                                                  
    sets_ <- paste(lapply(1:length(vars_),
                        function(x) 
                          paste("set ", vars_[x], " ?",cnt+x," ", sep="")
                  ), collapse="")    
    end_ <- "] ])"
    merged_ = paste(prev_, inp_, ask_, sets_, end_, sep="", collapse="")
    NLCommand(merged_, nl.obj=nl.obj)    
  }
  else if (is.vector(input)) {
    if (length(input) != NLReport(paste("count ",agentset,sep=" "), nl.obj=nl.obj)) {
      stop("Length of agentset not equal to length of input.")
    }
    if (length(var.name) != 1) {
      stop("For vector input you have to submit one agent variable name with argument var.name")
    }
    # construct processing string    
    start_ <- paste("(foreach sort ", agentset, sep="")
    inp_ <- paste("[",paste(input, collapse=" "), "]", collapse=" ")
    ask_ <- paste(" [ask ?1 [set ",var.name," ?2]])", sep="")
    merged_ <- paste(start_, inp_, ask_, sep=" ")
    NLCommand(merged_, nl.obj=nl.obj)    
  }
  else {
    stop("Input has to be a data.frame or vector.")  
  }
}


