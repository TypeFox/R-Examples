#' @export
.onAttach <- function(lib,pkg)
{
  packageStartupMessage("WARNING : Since Rancluster 0.92, the ranks have to be given to the package in the ranking notation (see convertRank function), with the following convention :
- missing positions are replaced by 0
- tied are replaced by the lowest position they share")
}
