`degreesOrKnots.fnc` <-
function(name) {
  s = strsplit(name, " ")[[1]][2]
  s2 = strsplit(s, "[,)]")
  return(as.numeric(s2[[1]][1]))
}

