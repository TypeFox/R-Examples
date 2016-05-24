ps.balance <- function(object,
                       sel           = NULL,
                       treat         = NULL,
                       stratum.index = NULL,
                       match.index   = NULL,
                       method        = "classical",
                       cat.levels    = 2,
                       alpha         = 5,
                       equal         = TRUE) 
{
  if (missing(object))
    stop("Argument 'object' is not given.")


  if (any(class(object)=="stratified.data.frame") |
      any(class(object)=="stratified.pscore") |
      any(class(object)=="matched.data.frame") |
      any(class(object)=="matched.data.frames") |
      any(class(object)=="matched.pscore") |
      any(class(object)=="data.frame"))
    
    UseMethod("ps.balance")
  else
    stop("Class of argument 'object' will not supported.")
 
}
