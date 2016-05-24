escape.gt.lt <- function(s) {
  s = gsub("<", ".$GElt$.", s, fixed=TRUE)
  s = gsub(">", ".$GEgt$.", s, fixed=TRUE)
  return (s)  
}
unescape.gt.lt <- function(s) {
  s = gsub(".$GElt$.", "<", s, fixed=TRUE)
  s = gsub(".$GEgt$.", ">", s, fixed=TRUE)
  return (s)  
}

trim_brackets <- function (x) gsub("^<+|>+$", "", x)

trim_space <- function (x) gsub("^\\s+|\\s+$", "", x)
