is.pemt <-
function(object) {
  if(!is(object, "pemt")) return(FALSE)
  len <- length(object)
  nomi <- sapply(object, names)
  if(!all(unlist(sapply(nomi[1:(len - 1)], "%in%", c("Tprobs", "Eprobs", "X", "Y"))))) return(FALSE)
  if(!all(nomi[[len]] %in% c("nk", "nomi", "coordsnames", "mpoints", "which.dire"))) return(FALSE)
  if (!all(sapply(object[1:(len - 1)], function(xx) is.multi_transiogram(xx$Tprobs) | is.null(xx$Tprobs)))) return(FALSE)
  return(TRUE)
}

