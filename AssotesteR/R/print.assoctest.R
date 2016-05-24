print.assoctest <-
function(x, ...)
{
  cat("\n", "\t", x$name, "\n\n")
  cat("Info:", "\n")
  print(x$args, quote=FALSE, print.gap=2)
  cat("\n")
  if (length(x) == 5) {
    if ("signs" %in% names(x) || "set" %in% names(x)) {
      print(round(unlist(x[1:2]), 6), print.gap=3)
    } else {
      print(round(unlist(x[1:3]), 6), print.gap=3)
    }
  } else print(round(unlist(x[1:2]), 6), print.gap=3)
  invisible(x)
}


