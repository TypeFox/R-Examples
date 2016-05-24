check_ped <- function(pedigree)
{
  if (ncol(pedigree) != 4) 
          stop("Pedigree must have 4 columns (id, m, f, obs).")

  # Check the format of pedigree - if all character may need to convert
  # to numeric. 
  pednum <- convertped(pedigree)

  nfounders <- sum(pednum[,2]==0 & pednum[,3]==0)
  if (nfounders !=4 & nfounders !=8) 
	stop("Cannot process a number of founders which is not 4 or 8")

  # Check that parents come before children. 
  sorted <- (pednum[,1]>pednum[,2] & pednum[,1]>pednum[,3])
  if (!any(sorted))
	stop("Pedigree is not sorted, parents should appear before children")

  return(pednum)
}

