NLGetAgentSet <-
#function(agent.var, agentset, as.data.frame=FALSE, df.col.names=NULL, nl.obj=NULL)
function(agent.var, agentset, as.data.frame=TRUE, agents.by.row=FALSE, as.vector=FALSE, nl.obj=NULL)
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
  
  # check for empty agentset
  #if (NLReport(paste("count",agentset),nl.obj=nl.obj) == 0) {    
  if (!(NLReport(paste("any? ",agentset),nl.obj=nl.obj))) {
    stop("The requested agentset is empty")
  }

  # create a vector
  if (as.vector == TRUE) {
    if (length(agent.var) != 1) {
      stop("as.vector=TRUE makes only sense if you request just one agent variable.")
    }
    avar <- c("map [[",agent.var,"] of ?] sort ", agentset)
    avar <- paste(avar, collapse="")
    resobj <- NLReport(avar, nl.obj=nl.obj)  
  }
  else {
    # create a data.frame
    if (as.data.frame == TRUE) { 
      str <- lapply(agent.var, function(x) {paste("NLReport(\"map [[",x,"] of ?] sort ",agentset,"\",nl.obj=nl.obj)",sep="")})
      str <- paste(str, collapse=",")
      str <- paste("resobj <- data.frame(",str,")",sep=" ")
      eval(parse(text=str))  
      names(resobj) <- agent.var
    }
    else {
      # create an "old-style" list
      if (agents.by.row == TRUE) {
        avar <- lapply(agent.var, function(x) {paste(c("[",x,"] of ?"), collapse="")} )
        avar <- c("map [(list ",avar,")] sort ", agentset)
      } 
      # create a "new-style" list
      else {
        avar <- lapply(agent.var, function(x) {paste(c("map [[",x,"] of ?] sort", agentset), collapse=" ")})
        avar <- paste(c("(list",avar,")"), collapse=" ")
      }
      avar <- paste(avar, collapse="")
      resobj <- NLReport(avar, nl.obj=nl.obj)
      if (agents.by.row == FALSE) {
        names(resobj) <- agent.var
      }
    }
  }  
  return (resobj)
}

