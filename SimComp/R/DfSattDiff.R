DfSattDiff <-
function(n,sd,type="Dunnett",base=1,ContrastMat=NULL) {


if (length(n)!=length(sd)) {
  stop("Lengths of n and sd must be equal")
}

if (!is.null(ContrastMat)) {
  if (ncol(ContrastMat)!=length(n)) {
    stop("Number of columns of ContrastMat and length of n must be equal")
  }
  C0 <- apply(X=ContrastMat, MARGIN=1, function(x) {
    all(x==0)
  })
  if (any(C0)) {
    cat("Warning: At least one row of ContrastMat is a vector with all components", "\n",
        "equal to zero", "\n")
  }
  Cmat <- ContrastMat
  type <- "User defined"
  if (is.null(rownames(Cmat))) {
    rownames(Cmat) <- paste("C", 1:nrow(Cmat), sep="")
  }
} else {
  type <- match.arg(type, choices=c("Dunnett", "Tukey", "Sequen", "AVE", "GrandMean", "Changepoint", 
    "Marcus", "McDermott", "Williams", "UmbrellaWilliams"))
  Cmat <- contrMat(n=n, type=type, base=base)
}
comp.names <- rownames(Cmat)

defrvec <- numeric(nrow(Cmat))
for (z in 1:nrow(Cmat)) {
defrvec[z] <- ( (sum((Cmat[z,])^2*sd^2/n))^2 ) / 
              sum( ( (Cmat[z,])^4*sd^4 ) / ( n^2*(n-1) ) ) }
names(defrvec) <- comp.names

return(defrvec)


}
