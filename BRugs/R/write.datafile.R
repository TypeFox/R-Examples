"write.datafile" <-
function (datalist, towhere, fill = TRUE){
  if (!is.list(datalist) || is.data.frame(datalist)) 
      stop("First argument to write.datafile must be a list.")
  cat(formatdata(datalist), file = towhere, fill = fill)
}
