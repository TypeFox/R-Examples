getHTTPErrorClass =
function(status)
{

  if("status" %in% names(status))
    status = status[["status"]]


  i = match(as.character(status), names(httpErrorClasses))

  if(is.na(i))
     "GenericHTTPError"
  else
     c(httpErrorClasses[i], "HTTPError")
}
