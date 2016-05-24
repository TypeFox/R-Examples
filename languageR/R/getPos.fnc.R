`getPos.fnc` <-
function(vec, pos) {
  if (pos == "end") return(length(vec))
  else {
    if (pos=="beg") return(1)
    else return(floor(length(vec)/2))
  }
}

