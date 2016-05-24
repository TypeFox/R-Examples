qordinal = function(p,mu,muAll)
## QORDINAL - finds the level of an ordinal variable from its quantile and a matrix of cumulative probabilities for j-1.
## Input arguments are:
## P, a vector of cumultative probabilities
## MU, an irrelevant list with cumulative probs of answer ($ETA1) and its previous value ($ETA2) as returned by predict("cumprob") function
## MUALL, a matrix of cumulative probabilities for j-1 across all levels j=1,...,J.
{
  rank=apply(muAll<=p,1,sum)
  levels = dimnames(muAll)[[2]] #levels are stored as column labels of muAll
  y = factor( levels[rank], levels=levels)
  return(y)
}