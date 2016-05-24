


model_comp <- function(object,...)
{
# object:  of class nrm or nlm
# ... : designs
  
  TYPE <- ifelse(class(object) == "nrm","NRM","NLM")
  
  dem <- list(...)
  
  #### COTROL & Q-MATRIX creation
  
  Qs <- lapply(dem,function(ed){

        design <- ctrl_design(design=ed,aDD=object$reshOBJ$aDD,gr=object$reshOBJ$gr,TYPE=TYPE)  
    
    # create Q matrix
    if(attr(object$reshOBJ$Qmat,"paraM") == "bock")
        {
          gdema <- grDMb(object$reshOBJ$aDD,object$reshOBJ$gr,design,TYPE=TYPE)
        } else if(attr(object$reshOBJ$Qmat,"paraM") == "01")
          {
          gdema <- grDM(object$reshOBJ$aDD,object$reshOBJ$gr,design,TYPE=TYPE)
          }  
    
   return(gdema) 
  })

  reshEX <- object$reshOBJ
  old <- reshEX$Qmat %*% object$etapar 
  
  if(class(object) == "nrm")
  {
  modres <- lapply(Qs,function(ds)
      {
        reshEX$Qmat <- ds # add new Q matrix
        newsv <- solve(t(ds) %*% ds) %*%  (t(ds) %*% old) # new starting values
        try(nrm(reshEX,etastart=newsv,ctrl=object$ctrl))
      })
  
    
  } else {
    
    
  modres <- lapply(Qs,function(ds)
  {
    reshEX$Qmat <- ds # add new Q matrix
    newsv <- solve(t(ds) %*% ds) %*%  (t(ds) %*% old) # new starting values
    try(nelm(reshEX,etastart=newsv,ctrl=object$ctrl))
  })
    
  }
  
  
  bigres <- list(object=object,modres=modres)
  class(bigres) <- "modc"
  


return(bigres)  
}  
  