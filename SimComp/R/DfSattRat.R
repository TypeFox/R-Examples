DfSattRat <-
function(n,sd,type="Dunnett",base=1,Num.Contrast=NULL,Den.Contrast=NULL,
                      Margin=NULL) {


if (length(n)!=length(sd)) {
  stop("Lengths of n and sd must be equal")
}

if ( is.null(Num.Contrast)+is.null(Den.Contrast)<2 ) {              # if at most one matrix is missing
  if (!is.null(Num.Contrast) & is.null(Den.Contrast)) {
    stop("Num.Contrast is specified, but Den.Contrast is missing")
  }
  if (is.null(Num.Contrast) & !is.null(Den.Contrast)) {
    stop("Den.Contrast is specified, but Num.Contrast is missing")
  }
  if (!is.null(Num.Contrast) & !is.null(Den.Contrast)) {
    if (nrow(Num.Contrast)!=nrow(Den.Contrast)) {
      stop("Number of rows of Num.Contrast and Den.Contrast must be equal")
    }
    if (ncol(Num.Contrast)!=length(n) | ncol(Den.Contrast)!=length(n)) {
      stop("Number of columns of Num.Contrast and Den.Contrast and number of groups must be equal")
    }
    NC0 <- apply(X=Num.Contrast, MARGIN=1, function(x) {
      all(x==0)
    })
    DC0 <- apply(X=Den.Contrast, MARGIN=1, function(x) {
      all(x==0)
    })
    if (any(c(NC0, DC0))) {
      cat("Warning: At least one row of Num.Contrast or Den.Contrast is a vector with all components", "\n",
          "equal to zero", "\n")
    }
    Num.Cmat <- Num.Contrast
    Den.Cmat <- Den.Contrast
    type <- "User defined"
    if (is.null(rownames(Num.Cmat)) && is.null(rownames(Den.Cmat))) {
      comp.names <- paste("C", 1:nrow(Num.Cmat), sep="")
    } else {
      if (any(rownames(Num.Cmat)!=rownames(Den.Cmat))) {
        comp.names <- paste(rownames(Num.Cmat), rownames(Den.Cmat), 
                       sep="/")
      } else {
        comp.names <- rownames(Num.Cmat)
      }
    }
  }
} else {                                                            # if both matrices are missing
  type <- match.arg(type, choices=c("Dunnett", "Tukey", "Sequen", "AVE", "GrandMean", "Changepoint", 
    "Marcus", "McDermott", "Williams", "UmbrellaWilliams"))
  Cmat <- contrMatRatio(n=n, type=type, base=base)
  Num.Cmat <- Cmat$numC
  Den.Cmat <- Cmat$denC
  comp.names <- Cmat$rnames
}

if (is.null(Margin)) {
  Margin <- rep(1,nrow(Num.Cmat))
}
if (is.numeric(Margin)) {
  if (length(Margin)==1) {
    Margin <- rep(Margin,nrow(Num.Cmat))
  }
  if (length(Margin)!=nrow(Num.Cmat)) {
    stop("Margin must be a single numeric value, or a numeric vector with length equal to the number of contrasts")
  }
}

defrvec <- numeric(nrow(Num.Cmat))
for (z in 1:nrow(Num.Cmat)) {
defrvec[z] <- ( (sum((Num.Cmat[z,]-Margin[z]*Den.Cmat[z,])^2*sd^2/n))^2 ) / 
              sum( ( (Num.Cmat[z,]-Margin[z]*Den.Cmat[z,])^4*sd^4 ) / ( n^2*(n-1) ) ) }
names(defrvec) <- comp.names

return(defrvec)


}
