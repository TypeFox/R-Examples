maximalvectors <- function(performances) {
  mvc <- mvc.create(performances)
  inds <- mvc.computeBESTindices(mvc)
  performances[inds,]
}

maximalvectors.indices <- function(performances) {
  mvc <- mvc.create(performances)
  return(mvc.computeBESTindices(mvc))
}

mvc.create <- function(perfMat) {
  mvc <- .jnew("fi/smaa/libror/r/MaximalVectorComputationRFacade",
                 as.vector(perfMat), as.integer(nrow(perfMat)))
  return(mvc)
}

mvc.computeBESTindices <- function(mvc) {
  mvind <- .jcall(mvc,
         "[I",
         method="computeBESTindices",
         simplify=TRUE)
  return(mvind)
}
