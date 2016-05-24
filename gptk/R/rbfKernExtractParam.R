rbfKernExtractParam <-
function (kern, only.values=TRUE,
                                 untransformed.values=TRUE) {
  params <- c(kern$inverseWidth, kern$variance)

  if ( !only.values )
    names(params) <- c("inverseWidth", "variance")

  return (params)
}
