#<<BEGIN>>
outm <- function(x,value="each",which.node=1)
#TITLE Output of Nodes
#DESCRIPTION
# Changes the output of Nodes
#KEYWORDS misc
#INPUTS
#{x}<<A \samp{mcnode} or a \samp{mc} object.>>
#[INPUTS]
#{value}<<The output of the \samp{mcnode} for multivariates nodes. May be "each" (default)
#if output should be provided for each variates considered independently, "none" for no output
#or a vector of name of function(s) (as a character string) that will be applied on the variates dimension
#before any output (ex: \samp{"mean"}, \samp{"median"}, \samp{c("min","max")}). The function should have no other arguments
#and send one value per vector of values (ex. do not use \samp{"range"}).>>
#{which.node}<<which node should be changed in a \samp{mc} object>>
#EXAMPLE
#data(total)
#total$xVUM2
### since outm = NULL
#summary(total$xVUM2)  
#x <- outm(total$xVUM2,c("min"))
#summary(x)

#VALUE
#\samp{x} with a modified \samp{outm} attribute.
#--------------------------------------------
{

  if(is.character(value)  && (value == "none" || value == "each" | all(sapply(value,exists,mode="function")))){

  if(inherits(x,"mcnode"))  attr(x,which="outm") <- value
      else if(inherits(x,"mc")) {
        if(is.character(which.node)) which.node <- which(names(x) %in% which.node)
        attr(x[which.node],which="outm") <- value}
    }
  else stop("value should be 'none','each' or the name a valid function")
  return(x)
}
