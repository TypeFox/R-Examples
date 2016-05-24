ssmrob <-
function(outcome, selection, control=heckitrob.control())
{
  type <- 0
  if (!inherits(selection, "formula")) {
    stop("argument 'selection' must be a formula")
  }
  if (length(selection) != 3) {
    stop("argument 'selection' must be a 2-sided formula")
  }
  if (inherits(outcome, "formula")) {
    if (length(outcome) != 3) {
      stop("argument 'outcome' must be a 2-sided formula")
    }
    type <- 2
  }
  else if (inherits(outcome, "list")) {
    if (length(outcome) == 1) {
      outcome <- outcome[[1]]
      type <- 2
    }
    else if (length(outcome) == 2) {
      if (inherits(outcome[[1]], "formula")) {
        if (length(outcome[[1]]) != 3) {
          stop("argument 'outcome[[1]]' must be a 2-sided formula")
        }
      }
      else stop("argument 'outcome[[1]]' must be either a formula or a list of two formulas")
      if (inherits(outcome[[2]], "formula")) {
        if (length(outcome[[2]]) != 3) {
          stop("argument 'outcome[[2]]' must be a 2-sided formula")
        }
      }
      else stop("argument 'outcome[[2]]' must be either a formula or a list of two formulas")
      type <- 5
    }
    else stop("argument 'outcome' must contain 1 or 2 components")
  }
  if(type==2)
  {
    return(heckitrob(outcome, selection, control=control))
  }
  else if(type==5)
  {
    return(heckit5rob(outcome[[1]], outcome[[2]], selection, control=control))
  }
}
