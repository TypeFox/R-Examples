format.TE <- function(TE, na=FALSE){
  TE <- meta:::rmSpace(TE)
  if (na) res <- format(TE)
  else res <- ifelse(is.na(TE), "", format(TE))
  res
}
