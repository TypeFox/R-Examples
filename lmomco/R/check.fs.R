"check.fs" <-
function(fs) {
  if(any(fs < 0) || any(fs > 1)) {
     warning("invalid nonexceedance probability")
     return(FALSE)
  }
  return(TRUE)
}

