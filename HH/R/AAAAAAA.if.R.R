"if.R" <-
function(r, s) {
  if (exists("is.R") && is.function(is.R) && is.R()) r
  else s
}

