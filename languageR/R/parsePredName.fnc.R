`parsePredName.fnc` <-
function(name) {
  s = strsplit(name, "[\\(\\)]")[[1]][2]
  s2 = strsplit(s, ", ")[[1]]
  return(list(baseName=s2[1], knots=as.numeric(s2[2])))
}

