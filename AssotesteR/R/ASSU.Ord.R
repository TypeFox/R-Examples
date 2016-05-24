ASSU.Ord <-
function(y, X, perm = 100)
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
  stat.assu = my_assu_method(U, V)
  assu.stat2 = stat.assu[3]  # stat ordered
  p2.assu = stat.assu[4]     # pval ordered
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)	
    ymean = mean(y)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      # center phenotype y
      y.perm = y[perm.sample] - ymean
      # get score vector
      U.perm = colSums(y.perm * X, na.rm=TRUE)
      perm.assu = my_assu_method(U.perm, V)
      x.perm[i] = perm.assu[4]		
    }
    # p-value 
    perm.pval = sum(x.perm < p2.assu) / perm  # ordered
  }
  
  ## results
  name = "ASSU.Ord: Adaptive Sum of Squared Score (Ordered)"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(assu.stat = assu.stat2,
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"
  return(res)
}

