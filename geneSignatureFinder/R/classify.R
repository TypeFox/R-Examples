classify <-
function(ddata) {
  if(is.matrix(ddata)) {
      ans <- pam(ddata, 2)
      ans$clusters <- ans$clustering
    } else {ans <- pamUnbiased(ddata)}
  return(ans)
}
