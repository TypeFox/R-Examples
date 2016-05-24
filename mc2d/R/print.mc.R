#<<BEGIN>>
print.mc <- function(x, digits=3,...)
#TITLE Prints a mcnode or a mc Object
#DESCRIPTION
# Print a description of the structure of the \samp{mc} or the \samp{mcnode} object.
#KEYWORDS print
#INPUTS
#{x}<<a \samp{mcnode} or a \samp{mc} object.>>
#[INPUTS]
#{digits}<<Number of digits to be used.>>
#{\dots}<<Further arguments to be passed to the print function.>>
#VALUE
#An invisible data frame.
#DETAILS
#SEE ALSO
# \code{\link{mcnode}} for \samp{mcnode} objects.
# \code{\link{mc}} for \samp{mc} objects.
#EXAMPLE

#CREATED 08-01-25
#REVISED 08-01-25
#--------------------------------------------
{
  outm <- lapply(x,attr,which="outm")
  type <- lapply(x,attr,which="type")
  dimm <- lapply(x,dim)
  nom <- names(x)

  sortie <- function(obj,x=1,out,typ){
    dimm <- dim(obj)
    rangem <- range(obj,na.rm=TRUE)
    return(data.frame(
      node=nom[i],
      mode=mode(obj),                            
      nsv=dimm[1],
      nsu=dimm[2],
      nva=nvariates,
      variate=x,
      min=rangem[1],
      mean=mean(obj,na.rm=TRUE),
      median=median(as.vector(obj) * 1,na.rm=TRUE),   # For logical
      max=rangem[2],
      Nas=sum(is.na(unclass(obj))),
      type=type[[i]],
      outm=out
    ))
  }
  
  res <- data.frame(NULL)
  for(i in seq_along(x)){
    nvariates <- dimm[[i]][3]
      if(outm[[i]][1]=="none") next
      if(outm[[i]][1]=="each"){
        for(j in 1:nvariates)
          res <- rbind(res,sortie(x[[i]][,,j,drop=FALSE],j,outm[[i]]))
      } else {
      for(func in outm[[i]]){
          tmp <- apply(x[[i]],c(1,2),get(func,mode="function"))
          res <- rbind(res,sortie(tmp,NA,func))
          }
  }
  }

 if(is.null(unlist(res))) {
  warning("Try to print a multivariate node with outm = none", call.=FALSE)
  res <- NULL}
  
 invisible(print(res,digits=digits,...))
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
print.mcnode <- function(x,...)
#ISALIAS print.mc
#--------------------------------------------
{
  print.mc(list(x=x),...)}
