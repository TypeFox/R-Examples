# methods for objects of class "cm" and "cmList"

summary.cmList <- function(object, ...)
{
  formulas <- lapply(lapply(object$results, formula), getx)
  ics <- lapply(object$results, IC.cm)
  res <- list()
  res$table <- cbind(IC=ics,
                     formula=formulas)
  res$ictype = object$results[[1]]$ictype
  colnames(res$table)[1] = res$ictype
  res$y = gety(formula(object$results[[1]]))
  
  class(res) <- "summary.cmList"
  return(res)
}

summaryTable <- function(object, index, ...)
{
  # TODO what should this summary table really look like for lme case?
  fit = getNthModel(object, index)
  do.lme = FALSE
  if (object$modeltype == "lme") {
    do.lme = TRUE
    tbl = round(summary(fit)$tTable, 5)
  } else {
    tbl = round(coef(summary(fit) ), 5)
  }
  tbl11 = tbl[1, 1]
  if (do.lme) {
    tbl.pct = cbind(round(tbl[-1, 1] / tbl11 * 100, 2), round(tbl[-1, 5], 3) )
  } else {
    tbl.pct = cbind(round(tbl[-1, 1] / tbl11 * 100, 2), round(tbl[-1, 4], 3) )
  }
  colnames(tbl.pct) = c("Estimate (%)", "p-Value")
  tbl.pct = rbind(c(tbl11, 0), tbl.pct)
  return(tbl.pct)
}

getNthModel <- function(object, index)
{
  if(index == 1) {
    return(getFirstModel(object))
  }
  # otherwise, need to calculate the model first
  tmp.gds=object$data
  random=object$random
  modeltype=object$modeltype
  return(with(tmp.gds,eval(object$results[[index]]$call)))
}

getFirstModel <- function(object)
{
  return(object$first$fit)
}

gety <- function(formula)
{
  return(deparse(formula[[2]], width.cutoff=500))
}
getx <- function(formula)
{
  return(deparse(formula[[3]], width.cutoff=500))
}

print.summary.cmList <- function(x, ...)
{
  cat(paste("models sorted by ", x$ictype[[1]], ", first model is smallest:\n"))
  cat("y: \"", paste(x$y[[1]],sep=""),"~\"\n")
  print(x$table, print.gap=2, digits=2, ...)
}

IC.cm <- function(object)
{
  return(object$IC)
}

formula.cm <- function(x, ...)
{
  return(x$formula)
}

summary.cm <- function(object, ...)
{
  formula <- deparse(object$formula, width.cutoff=500)
  ic <- list(object$IC)
  res <- list()
  res$table <- cbind(IC=ic,
                     formula=formula)
  colnames(res$table)[1] = object$ictype
  class(res) <- "summary.cm"
  return(res)
}

print.summary.cm <- function(x, ...)
{
  print(x$table, print.gap=2, digits=2, ...)
}
