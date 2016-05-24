# power for detecting interaction effect in 2-way ANOVA
# nTotal = number of observations in total
# a = number of levels in factor 1
# b = number of levels in factor 2
# effsize = effect size
# alpha = type I error rate
# nTests = number of tests if multiple testing

powerInteract=function(nTotal, a, b, effsize, alpha=0.05, nTests=1)
{
  alpha2=alpha/nTests

  nPerCell=floor(nTotal/(a*b))
  df1=(a-1)*(b-1)
  df2=a*b*(nPerCell-1)

  F0=qf(p=1-alpha2, df1=df1, df2=df2, ncp=0)

  ncp= nPerCell*effsize^2
  power=1-pf(q=F0, df1=df1, df2=df2, ncp=ncp)
  return(power)
}

