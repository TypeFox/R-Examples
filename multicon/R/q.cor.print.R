q.cor.print <-
function(obj, var.content, initial, rnd=2, EXPORT=FALSE, short=FALSE) {
  .Deprecated("print.q.cor", msg = "q.cor.print is deprecated. Please use the print function.")
  answer <- print(x=obj, var.content = var.content, initial = initial, rnd = rnd, EXPORT = EXPORT, short = short)
  return(answer)
}
