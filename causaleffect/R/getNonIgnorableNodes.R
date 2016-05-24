getNonIgnorableNodes <-
function(x) {
  var <- c()
  if (x$recursive) {
    for (i in 1:length(x$children)) var <- c(var, getNonIgnorableNodes(x$children[[i]]))
  } else {
    var <- c(var, x$var, x$cond)
  }
  return(var)
}
