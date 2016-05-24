additional_comp <-
function(by, svydat)  
{

  
  #########################################################
  ###########        N cases       ######################
  #########################################################   
  
  tabnamsplit <- all.vars(by)
  aggcom      <- paste("list(",paste0("svydat$variables$",tabnamsplit,collapse=","),")",sep="")
  
  Ncases      <- aggregate(svydat$variables[,1], eval(parse(text=aggcom)), FUN=length)
  
  #########################################################
  ###########     Sum of weights    ####################
  #########################################################   
  
  Sumweights <- aggregate(svydat$pweights, eval(parse(text=aggcom)), FUN=sum)

  return(list(Ncases=Ncases, Sumweights=Sumweights))
  
}
