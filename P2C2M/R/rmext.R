rmext <-
function (instring, delimiter="\\.") {
  # function for returning instring without extension (i.e. ".txt" or ".trees")
  # get all instances of "."
  alist = stringr::str_locate_all(instring, delimiter)
  # get position of last element
  pos2 = tail(alist[[1]], 1)[1]-1
  return(substr(instring, 1, pos2))
}
