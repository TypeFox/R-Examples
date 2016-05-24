groupSSANearestNeighbour = function(
##title<< Group SSA eigentriples by finding nearest euclidian neighbours
  x          ##<< object of class ssa (e.g. the results from a call to ssa)
  , ...      ##<< other objects that can be passed to the function but which
             ##   are not used. This is only implemented to make the function
             ##    identical in its call to the grouping.auto function.
)
##description<<
##This function finds groups in SSA eigentriples by reconstructing
##individual eigentriples back to their original (time)  domain and by
##determining the nearest euclidian neighbour for each eigentriple.
##Groups with more than two members are constructed by identifying
##groups with a very similar Euclidian distance.
##seealso<<
##\code{\link{ssa}}, \code{\link{grouping.auto}}
{
  r.grouping       <- reconstruct(x)
  n.comp           <- length(r.grouping)
  r.matrix.group   <- matrix(unlist(r.grouping[1:n.comp]), nrow = n.comp, byrow=TRUE)
  dist.matrix      <- as.matrix(dist(r.matrix.group, upper = TRUE))
  parts            <- 1:n.comp
  groups.ssa       <- list()
  step             <- 1
  for (i in 1:n.comp)
  {
    part1       <- parts[1]
    parts       <- setdiff(parts, part1)
    part2       <- parts[which(dist.matrix[parts, part1] < 0.8 * mean(dist.matrix[parts, part1]))]
    groups.ssa[[step]] <- c(part1, part2)
    parts       <- setdiff(parts, part2)
    step        <- step + 1
    if (length(parts) == 0)
      break
  }
  ##value<< list: list indicating the grouping of the SSA eigentriples
  groups.ssa
}
