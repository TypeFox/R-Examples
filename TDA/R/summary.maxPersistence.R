summary.maxPersistence <-
function(object, ...) {
  call <- attributes(object)[["call"]]
  maxNum <- max(object[["sigNumber"]])
  maxNumParam <- object[["parameters"]][
      which(object[["sigNumber"]] == max(object[["sigNumber"]]))]
  maxSig <- max(object[["sigPersistence"]])
  maxSigParam <- object[["parameters"]][
      which(object[["sigPersistence"]] == max(object[["sigPersistence"]]))]
  out <- list(call = call, maxNum = maxNum, maxNumParam = maxNumParam,
              maxSig = maxSig, maxSigParam = maxSigParam)
  class(out) <- c("summary.maxPersistence")
  return(out)
}