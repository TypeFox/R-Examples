ps.estimate <- function(object,
                        resp,
                        treat         = NULL,
                        stratum.index = NULL,
                        match.index   = NULL,
                        adj           = NULL,
                        weights       = "rr",
                        family        = "gaussian",
                        regr          = NULL,
                        ...
                        ) 
{
  if (missing(object))
    stop("Argument 'object' is not given.")


  if (any(class(object)=="stratified.data.frame") |
      any(class(object)=="stratified.pscore") |
      any(class(object)=="matched.data.frame") |
      any(class(object)=="matched.data.frames") |
      any(class(object)=="matched.pscore") |
      any(class(object)=="data.frame"))
    
    UseMethod("ps.estimate")
  else
    stop("Class of argument 'object' will not supported.")
 
}
