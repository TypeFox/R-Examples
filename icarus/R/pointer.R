# copyright (C) 2015 A.Rebecq

newPointer <- function(inputValue) 
  { 
  object=new.env(parent=globalenv()) 
  object$value=inputValue 
  class(object)='pointer'
  
  return(object) 
}  