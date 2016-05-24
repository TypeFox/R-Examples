print.bart_spher <- function(x, ...){
  if(interactive()) writeLines("")
  writeLines("\tBartlett's Test of Sphericity\n")
  writeLines(paste0("Call: ", deparse(x$call), "\n"))
  colns <- format(c("X2","df","p-value"), justify = "right")
  writeLines(paste0(colns[1L], " = ", round(x$X2, 3L)))
  writeLines(paste0(colns[2L], " = ", round(x$df, 0L)))
  pv <- my.format.pval(x$p.value, 5L)
  writeLines(paste0(colns[3L], ifelse(pure_all.equal(0, round(x$p.value, 5L)), " ", " = "), pv, collapse=""))
  if(interactive()) writeLines("")
  if(x$warn){
    warning(paste0("Used n = ", round(x$n,2L), "."), call. = FALSE, immediate. = TRUE)
  }
}
