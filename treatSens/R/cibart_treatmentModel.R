probitEM <- function(maxBackstepIterations = 30L) {
  if (maxBackstepIterations < 0) stop('illegal max iterations')
  structure(namedList(maxIter = maxBackstepIterations), class = "probitEMTreatmentModel")
}

bart <- function(k = 2, ntree = 50, keepevery = 10)
{
  if (k <= 0.0 || ntree <= 0 || keepevery <= 0)
    stop('illegal binary bart option')
  
  structure(namedList(k = as.double(k), ntree = as.integer(ntree), keepevery = as.integer(keepevery)),
            class = "bartTreatmentModel")
}

probitStudentTPrior <- function(df = 3, scale = 4.0) {
  if (df <= 0.0) stop('illegal prior degrees of freedom')
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(df, scale, family = "studentt"), class = "probitTreatmentModel")
}

probitCauchyPrior <- function(scale = 4.0) {
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(df = 1, scale, family = "studentt"), class = "probitTreatmentModel")
}

probitNormalPrior <- function(scale = 4.0) {
  if (scale <= 0.0) stop('illegal prior scale')

  structure(namedList(scale, family = "normal"), class = "probitTreatmentModel")
}

probit <- function(family = "cauchy", ...) {
  matchedCall <- match.call()
  familyMissing <- missing(family)
  
  if (!is.null(family) && !(family %in% c("normal", "flat", "cauchy", "t")))
    stop('unsupported family type')

  if (is.null(family) || identical(family, "flat")) return(structure(list(family = "flat"), class = "probitTreatmentModel"))
  
  familyCall <- matchedCall
  matchIndices <- match(names(familyCall), "family")

  if (!familyMissing)
    if (any(!is.na(matchIndices))) {
      familyCall <- familyCall[is.na(matchIndices)]
    } else {
      familyCall <- familyCall[-2]
    }
  
  familyCall[[1]] <-
    switch(family,
           cauchy = quoteInNamespace(probitCauchyPrior),
           t      = quoteInNamespace(probitStudentTPrior),
           normal = quoteInNamespace(probitNormalPrior))

  eval(familyCall, parent.frame(1))
}

