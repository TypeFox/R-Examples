assignList <- function(aList, pos = -1, envir = as.environment(pos), inherits = FALSE){
  if(is.null(nms <- names(aList))) stop("names(aList) is NULL")
  if(any(nms == "")) stop("blank name")
  if(missing(envir) && pos < 0)
    envir <- parent.frame(-pos)
  for(nm in nms)
    assign(x = nm, value = aList[[nm]], envir = envir, inherits = inherits)
}
