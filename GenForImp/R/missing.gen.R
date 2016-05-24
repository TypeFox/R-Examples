missing.gen <-
function(mat, nummiss){
  p <- ncol(mat)
  repeat{
  mmiss <- missing.gen0(mat, nummiss)
  ind.na <- is.na(mmiss)
  max.na <- max(as.numeric(names(table(apply(ind.na, 1, function(x) table(x)["TRUE"])) ))) 
  if(max.na < p) {break}
  }
  mmiss
}
