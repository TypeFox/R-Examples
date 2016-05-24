getLabel <- function (x, withnames = TRUE)
{
    cl <- as.character(class(x)[1])
    slots <-  slotNames(param(x))
    slots <-  slots[slots != "name"]
    qparamstring <-""
    nrvalues <-  length(slots)
    if(nrvalues > 0){
          values <-  numeric(nrvalues)
      for(i in 1:nrvalues)
        values[i] <-  attributes(attributes(x)$param)[[slots[i]]]
    if( withnames)
        nparamstring <-  paste(slots, "=", values, collapse = ", ")
   else
        nparamstring <-  paste(values, collapse = ", ")
    qparamstring <- paste("(",nparamstring,")",sep="")
  return(paste(cl,qparamstring,sep=""))
  }
} 
