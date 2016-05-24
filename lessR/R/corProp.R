corProp <- 
function (R=mycor, 
          main=NULL, heat.map=TRUE, bottom=3, right=3, 
          pdf.file=NULL, pdf.width=5, pdf.height=5) {


  # cor matrix:  mycor as class out_all, mycor$cors, or stand-alone matrix
  cor.nm <- deparse(substitute(R))
  .cor.exists(cor.nm)  # see if matrix exists in one of the 3 locations
  if (class(R) == "out_all")
    R <- eval(parse(text=paste(cor.nm, "$cors", sep="")))  # go to $cors 

  NVOld <- as.integer(nrow(R))

  Label <- integer(length=NVOld)
  NVC <- NVOld
  Diagon <- as.integer(0)
  Power <- as.integer(0)


  # Get the float version of Power, RPower
  RPower <- ifelse (Power == 0, 1.0, as.numeric(Power) / 100.0) 

  # Compute sum of squares for each column and store in Diag
  CP <- crossprod(R)
  Diag <- diag(CP)

  # Compute sum of cross-products and store in upper triangle.  Leave
  #    the correlations in the lower triangle with diagonal unchanged
  #    R[I,J) is the cross-product, R[J,I) is the original correlation

  for (I in 1:(NVC-1)) {
    for (J in (I+1):NVC) {
      R[I,J] <- 0
      for (K in 1:NVC) {
        if (I <= K) X1 <- R[K,I]
        if (I > K) X1 <- R[I,K]
        if (J <= K) X2 <- R[K,J]
        if (J > K) X2 <- R[J,K]
        R[I,J] <- R[I,J] + X1*X2
      }
    }
  }

  # Normalize cross products, i.e., obtain the proportionality coefs
  #   excluding the diagonal.  If the diagonal is ignored, Diagon<-0, then
  #   correspondingly reduce the appropriate sums of squares and
  #   cross-products

  for (I in 1:(NVC-1)) {
    for (J in (I+1):NVC) {
      if (Diagon == 0) {
        RII <- R[I,I]
        RJJ <- R[J,J]
        RJI <- R[J,I]
        D1 <- Diag[I] - (RII**2+RJI**2)
        D2 <- Diag[J] - (RJJ**2+RJI**2)
        R[I,J] <- R[I,J] - ((RII*RJI) + (RJJ*RJI))
      }
      else if (Diagon == 1) {
        D1 <- Diag(I)
        D2 <- Diag(J)
      }
      R[I,J] <- R[I,J] / (sqrt(D1*D2))
      R[I,J] <- R[I,J]**RPower
      R[J,I] <- R[I,J]
    }
  }

  # Set the diagonal to 1.OO
    for (I in 1:NVC) R[I,I] <- 1.00

  # assign names
  nm <- character(length=NVOld)
  nm <- dimnames(R)[[1]]
  dimnames(R) <- list(nm, nm)


  if (heat.map) {
    if (is.null(main)) main <- "Item Proportionalities"
   .corcolors(R, NVOld, main, bottom, right, diag=0,
              pdf.file, pdf.width, pdf.height)
  }

  # finish
  cat("\n")
  return(round(R,2))

}
