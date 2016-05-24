findWidgetByType <-
  #
  # Recursively search a widget tree for the first occurence of
  # a widget of the specified type.
  #
function(win, gType = "GtkMenuBar", verbose = FALSE)
{
 if(verbose)
   print(class(win))

 if(is.function(gType)) {
   if(gType(win))
     return(win)
 } else if(as.character(gType) %in% class(win)) {
     return(win)
 }

 if("GtkContainer" %in% class(win)) {
  for(i in win$GetChildren()) {
   tmp = findWidgetByType(i, gType, verbose = verbose)
   if(!is.null(tmp))
     return(tmp)
  }
 }

 return(NULL)
}
