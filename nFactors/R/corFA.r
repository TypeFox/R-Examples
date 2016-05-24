"corFA" <-
function(R, method="ginv") {
 R <- as.matrix(R)
 if (method == "ginv") return(R - ginv(diag(diag(ginv(R)))))
 }
