get.expression <-
function(x) {
  P <- ""
  s.print <- length(x$sumset) > 0
  sum.string <- ""
  do.string <- NULL
  cond.string <- ""
  if (x$fraction) P <- "\\frac{"
  if (s.print) {
    sum.string <- paste(x$sumset, sep = "", collapse = ",")
    P <- paste(P, "\\sum_{", sum.string, "}[", sep = "", collapse = "")
  }
  if (x$recursive) {
    for (i in 1:length(x$children)) P <- paste(P, get.expression(x$children[[i]]), sep = "", collapse = ",")
  } else {
    var.string <- paste(x$var, sep = "", collapse = ",")
    if (x$star) P <- paste(P, "P^*(", var.string, sep = "", collapse = "")
    else if (length(x$domain) > 0) P <- paste(P, "P^{(", x$domain, ")}(", var.string, sep = "", collapse = "")
      else P <- paste(P, "P(", var.string, sep = "", collapse = "")
    if (length(x$do) > 0) {
      do.string <- paste("do(", x$do, ")", sep = "")
    }
    if (length(x$cond) > 0) cond.string <- paste(x$cond, sep = "", collapse = ",")   
    if (length(x$cond) > 0 | length(x$do) > 0) {
      cond.string <- c(do.string, x$cond)
      cond.string <- paste(cond.string, sep = "", collapse = ",")
      cond.string <- paste("\\vert ", cond.string, ")", sep = "", collapse = "") 
    }
    else cond.string <- ")"
    P <- paste(P, cond.string, sep = "")
  }
  if (s.print) P <- paste(P, "]", sep = "", collapse = ",")
  if (x$fraction) { 
    P <- paste0(P, "}{")
    P <- paste(P, get.expression(x$divisor), sep = "", collapse = ",")
    P <- paste0(P, "}")
  }
  return(P)
}
