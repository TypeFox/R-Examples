


rv.any.na <- function (x) {
  # NAME
  #  rv.any.na - Which components have missing values?
  #
  if (is.rvsummary(x)) {
    return(unlist(rvattr(x, "NAS"))>0)
  } else {
    return(colSums(is.na(sims(as.rv(x))))>0)
  }
}



rv.all.na <- function (x)
{
  # NAME
  #  rv.all.na - Which components are completely missing?
  if (is.rvsummary(x)) {
    return(unlist(rvattr(x, "NAS"))==1)
  } else {
    return(colSums(is.na(sims(as.rv(x))))==1)
  }
}

