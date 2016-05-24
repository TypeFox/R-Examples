`getRange.fnc` <-
function(lst) {
  v = vector()
  for (i in 1:length(lst)) {
    if (is.data.frame(lst[[i]])) {
      if ("lower" %in% colnames(lst[[i]])) {
        v = c(v, as.vector(lst[[i]][,c("Y", "lower", "upper")]))
      } else {
        v = c(v, as.vector(lst[[i]][,"Y"]))
      }
    } else {
      for (j in 1:length(lst[[i]])) {
        if ("lower" %in% colnames(lst[[i]][[j]])) {
           v = c(v, as.vector(lst[[i]][[j]][,c("Y", "lower", "upper")]))
        } else {
           v = c(v, as.vector(lst[[i]][[j]][,"Y"]))
        }
      }
    }
  }
  return(range(v))
}

