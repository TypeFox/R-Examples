washout <-
function(x, washcol, method = "index") {
  if (method == "value") {
    washvec = cut(x,length(washcol),labels = FALSE)
  } else if (method == "index") {
    washvec = cut(1:length(x), length(washcol), labels = FALSE)
  } else {break}
  washout = washcol[washvec]
  washout
}
