## #########################################################
## Function to stratifiy data dependent of the class of data
## #########################################################

ps.makestrata <- function(object,
                          breaks             = NULL, 
                          name.stratum.index = "stratum.index",  
                          stratified.by      = NULL,
                          ...
                          )
{
  if (missing(object))
    stop("Argument 'object' is not given.")
  
  if (any(class(object)=="pscore") |
      any(class(object)=="data.frame"))
    
    UseMethod("ps.makestrata")
  else
    stop("Class of argument 'object' will not supported.")
  
}
