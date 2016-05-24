whiteKernExtractParam <-
function (kern, only.values=TRUE,
                                   untransformed.values=TRUE) {
  params <- c(kern$variance)

  if ( !only.values ) {
    names(params) <- c("variance")
  }

  return (params)
}
