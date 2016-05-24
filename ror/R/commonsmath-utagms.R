utagms.create <- function(perfMat) {
  model <- .jnew("fi/smaa/libror/r/UTAGMSSolverRFacade", as.vector(perfMat), as.integer(nrow(perfMat)))
  list(model=model,rownames=rownames(perfMat),colnames=colnames(perfMat))
}

utagms.solve <- function(ror) {
  val <- .jcall(ror$model, "I", method="solve")
  return(as.logical(val))
}

utagms.printModel <- function(ror, necessary, aind, bind) {
  aind = aind-1
  bind = bind-1
  .jcall(ror$model, "V", method="printModel", as.logical(necessary), as.integer(aind),
         as.integer(bind))
}

utagms.setStrictValueFunctions <- function(ror, strict) {
  .jcall(ror$model, "V", method="setStrictValueFunctions", as.logical(strict))
}

utagms.getNecessaryRelation <- function(ror) {
  rel <- .jcall(ror$model, "[[D", method="getNecessaryRelation", simplify=TRUE)
  rownames(rel) <- ror$rownames
  colnames(rel) <- ror$rownames
  return(rel)
}

utagms.getPossibleRelation <- function(ror) {
  rel <- .jcall(ror$model, "[[D", method="getPossibleRelation", simplify=TRUE)
  rownames(rel) <- ror$rownames
  colnames(rel) <- ror$rownames
  return(rel)
}

