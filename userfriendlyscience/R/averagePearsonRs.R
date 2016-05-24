averagePearsonRs <- function(rs, ns, FishersZ=TRUE) {
  if (FishersZ) {
    return(convert.fisherz.to.r(averageFishersZs(convert.r.to.fisherz(rs), ns)));
  } else {
    warning("Sorry, Alexander's method (1990, Bulletin of the Psychonomic Society) not yet implemented!");
  }
}