ASSUW.Ord <-
function(y, X, perm=100)
{
  ## checking arguments
  Xy_perm = my_check(y, X, perm)
  y = Xy_perm$y
  X = Xy_perm$X
  perm = Xy_perm$perm
  
  ## get U and V
  getuv = my_getUV(y, X)
  U = getuv$U
  V = getuv$V
  ## run score method
  stat.assuw = my_assuw_method(U, V)
  assuw.stat2 = stat.assuw[3]  # stat ordered
  p2.assuw = stat.assuw[4]     # pval ordered
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    p2.perm = rep(0, perm)
    ymean = mean(y)	
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      U.perm = colSums(y.perm * X, na.rm=TRUE)
      perm.assuw = my_assuw_method(U.perm, V)
      p2.perm[i] = perm.assuw[4]		
    }
    # p-value 
    perm.pval = sum(p2.perm < p2.assuw) / perm  # ordered
  }
  
  ## results
  name = "ASSUW.Ord: Adaptive Weighted Sum of Squared Score (Ordered)"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(assuw.stat = assuw.stat2, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

