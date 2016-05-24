print.ellipsemake <-
function(x,...) {
print(x$values[c("area","lax","retention","coercion")])
invisible(x)}
