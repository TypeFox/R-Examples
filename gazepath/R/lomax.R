lomax <-
function(x){
  which(diff(c(TRUE, diff(x) >= 0, FALSE)) < 0)
}
