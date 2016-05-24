`getKnots.fnc` <-
function(colnms, xlb) {
  pos = grep(paste("rcs\\(", xlb, ",", sep=""), colnms)
  tmp = strsplit(colnms[pos[1]], ")")[[1]][1]
  return(as.numeric(strsplit(tmp, " ")[[1]][2]))
}

