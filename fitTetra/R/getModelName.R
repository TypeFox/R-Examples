getModelName <-
function  (mutype,ptype) {
  s <- getMuModelName(mutype)
  if (ptype=="p.HW") {
    s <- paste(s,"HW")
  }
  s
}
