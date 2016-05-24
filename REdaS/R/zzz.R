my.format.pval <- function(p, digits){
  gsub("^0\\.", ".", format.pval(pv = round(p, digits), digits = digits, scientific = 7L))
}

pure_all.equal <- function(tgt, cur){
  return(isTRUE(all.equal.numeric(tgt, cur, check.attributes = FALSE)))
}
