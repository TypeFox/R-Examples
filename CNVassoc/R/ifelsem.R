ifelsem <-
function(test, yes, no){
  out <- no
  out[test, ] <- yes[test, ]
  out
}

