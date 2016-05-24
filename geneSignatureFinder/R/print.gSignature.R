print.gSignature <-
function(x, ...) {
  message(paste("Signature name:", x$signatureName))
  message(paste("starting gene(s): ", x$startingSignature, 
                ", with log(p.value, 10) = ", round(log(x$startingPValue, 10), 3), sep = ""))
  message(paste("signature of ", length(x$signature), " elements with log(p.value, 10) = ", round(log(x$pValue, 10), 3), sep = ""))
  message("\t", paste(x$signature, collapse = ", "))
#  NextMethod("print")
}
