`A.chc` <-
function(hc) {
  hts <- diff(c(0,hc$hei))
  unA <- .Call(UA,hc$mer,hts)
  return(unA/sum(hts))
}
