## #####################################################
## Function to match data dependent of the class of data
## #####################################################

ps.match <- function(object,
                     object.control     = NULL,
                     matched.by         = NULL,   
                     control.matched.by = matched.by,
                     who.treated        = 1,
                     treat              = NULL,
                     name.match.index   = "match.index",
                     ratio              = 1,
                     caliper            = "logit",
                     x                  = 0.2,
                     givenTmatchingC    = TRUE,
                     bestmatch.first    = TRUE,
                     setseed            = FALSE,
                     combine.output     = TRUE)
{
  if (missing(object))
    stop("Argument 'object' is not given.")
  
  if (any(class(object)=="pscore") |
      any(class(object)=="data.frame"))
    
    UseMethod("ps.match")
  else
    stop("Class of argument 'object' will not supported.")
}
