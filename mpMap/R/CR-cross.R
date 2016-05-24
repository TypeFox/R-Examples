CR_cross <-
function(mpcross)
{
  nfounders <- nrow(mpcross$founders)
  nfinals <- nrow(mpcross$finals)
  nped <- nrow(mpcross$pedigree)
  
  funnel <- .C("pf", as.integer(mpcross$id), as.integer(mpcross$pedigree[,1]), as.integer(mpcross$pedigree[,2]), as.integer(mpcross$pedigree[,3]),  as.integer(nfinals), as.integer(nfounders), as.integer(nped), out=integer(length=nfinals*nfounders), PACKAGE="mpMap")
  funnel$out <- match(funnel$out, mpcross$fid)

  return(matrix(data=funnel$out, nrow=nfinals, ncol=nfounders, byrow=TRUE))
}

