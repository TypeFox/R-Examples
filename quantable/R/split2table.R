#' splits names and creates a matrix
#' @export
#' @param names vector with names
#' @param split patter to use to split
#' @return matrix
#'
#' @examples
#' dat = c("bla_ra0/2_run0","bla_ra1/2_run0","bla_ra2/2_run0")
#' split2table(dat,split="\\_|\\/")
split2table <- function(names,split="\\||\\_")
{
  cnamessplit <- strsplit(names,split)
  protnam <- matrix(NA, ncol=length(cnamessplit[[1]]),nrow=length(cnamessplit))
  print(dim(protnam))
  for(i in 1:length(cnamessplit))
  {
    protnam[i,] <- cnamessplit[[i]]
  }
  return(protnam)
}
