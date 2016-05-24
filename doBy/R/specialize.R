specialize <- function(FUN, arglist){
  expr1 <- as.expression(body(FUN))
  expr2 <- do.call("substitute", list(expr1[[1]], arglist))
  gg  <- formals(FUN)
  idx <- match(names(arglist), names(gg))
  idx <- idx[!is.na(idx)]
  if (length(idx)>0){ gg  <- gg[-idx]}
  as.function(c(gg, expr2))
}
