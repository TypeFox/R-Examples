`getRoot.fnc` <-
function(xlabel) {
  a = regexpr("poly\\(", xlabel)
  if (a==1) {
    base = substr(xlabel, attr(a,"match.length"),nchar(xlabel))
    tmp = strsplit(base, ",")[[1]][1]
    return(substr(tmp, 2, nchar(tmp)))
  } else
    base=xlabel
  return(base)
}

