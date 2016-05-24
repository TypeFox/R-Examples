print.penalty <-
function (x, ...)
{
  cat (paste ("\nPenalty = ", x$penalty, " (lambda = ", x$lambda, ")\n", sep = ""))
  invisible (x)
}

