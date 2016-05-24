library(allelic)
pvalue <- allelic.exact.test(0,0,0,0,0,1)
if (pvalue != 1)
  stop("1: pvalue should be equal to 1");

pvalue <- allelic.exact.test(160,80,60,160,160,30)
r <- round(pvalue * 100000)
if ( r != 50743 )
  stop("2: bad pvalue, should be rounded to 0.50743")
