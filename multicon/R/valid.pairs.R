valid.pairs <-
function(x, y) {
  Tot <- nrow(data.frame(x,y))
  Miss <- sum(is.na(rowSums(data.frame(x,y))))
  Valid <- Tot - Miss
  Pct <- Valid / Tot
  out <- list("Tot"=Tot,"Miss"=Miss,"Valid"=Valid,"Pct"=Pct)
  return(out)
}
