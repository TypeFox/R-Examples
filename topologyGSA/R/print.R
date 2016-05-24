print.clique.var.test <- function(x, ...) {
  cat0("\n        Clique Variance Test\n\n")

  ns <- x$var.names
  cat0("data: ", ns[[1]], ", ", ns[[2]], " and ", ns[[3]], "\n\n")

  cat0("cliques:\n")
  for (i in seq_along(x$lambda.value)) {
    cat0(i, ": ",
         "lambda = ", x$lambda.value[i], ", ",
         "p-value = ", x$p.value[i], ", ",
         "equal variances: ", x$var.equal[i], "\n")
  }

  cat0("\n")
}

print.clique.mean.test <- function(x, ...) {
  cat0("\n        Clique Mean Test\n\n")

  ns <- x$var.names
  cat0("data: ", ns[[1]], ", ", ns[[2]], " and ", ns[[3]], "\n")
  cat0("paired: ", !("lambda.value" %in% names(x)), "\n\n")

  cat0("cliques:\n")
  for (i in seq_along(x$lambda.value)) {
    cat0(i, ": ",
         "t = ", x$t.value[i], ", ",
         "p.value = ", x$p.value[i],
         "\n")
  }

  cat0("\n")
}


print.pathway.var.test <- function(x, ...) {
  cat0("\n        Pathway Variance Test\n\n")

  ns <- x$var.names
  cat0("data: ", ns[[1]], ", ", ns[[2]], " and ", ns[[3]], "\n\n")

  cat0("lambda = ", x$lambda.value, ", ",
       "df = ", x$df, ", ",
       "p-value = ", x$p.value, ", ",
       "equal variances: ", x$var.equal, "\n\n")
}

print.pathway.mean.test <- function(x, ...) {
  cat0("\n        Pathway Mean Test\n\n")

  ns <- x$var.names
  paired <- !("lambda.value" %in% names(x))

  cat0("data:  ", ns[[1]], ", ", ns[[2]], " and ", ns[[3]], "\n")
  cat0("paired: ", paired, "\n\n")

  cat0("t = ", x$t.value, ", p-value = ", x$p.value)
  if (!paired)
    cat0(", ",
         "df.mean = ", x$df.mean, ", ",
         "equal variances: ", x$var.equal)

  cat0("\n\n")
}


cat0 <- function(...) cat(..., sep="")
