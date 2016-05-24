# findindexM function, called in selectterms

findindexM = function(indexM, nterms)
{
  coeff = 0
  termindex = length(indexM)
  i = 0
  j = 0
  k = 0
  l = 0
  m = indexM[termindex]
  
  if (termindex > 1)
    l = indexM[termindex-1]
  if (termindex > 2)
    k = indexM[termindex-2]
  if (termindex > 3)
    j = indexM[termindex-3]
  if (termindex > 4)
    i = indexM[termindex-4]
  
  for (p in 1:i-1)
    coeff <- coeff + choose(nterms-p, 4)
  print(coeff)
  for (p in 1:j-i-1)
    coeff <- coeff + choose(nterms-i-p, 3)
  print(coeff)
  for (p in 1:k-j-1)
    coeff <- coeff + choose(nterms-j-p, 2)
  print(coeff)
  for (p in 1:l-k-1)
    coeff <- coeff + choose(nterms-k-p, 1)
  print(coeff)
  
  coeff = coeff + m - l
  
  return(coeff)
}