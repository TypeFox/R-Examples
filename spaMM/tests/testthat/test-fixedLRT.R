cat("\ntest fixedLRT:")
# fixedLRT

data(blackcap)
## result comparable to the corrHLfit examples based on blackcap
fl <- fixedLRT(null.formula=migStatus ~ 1 + Matern(1|latitude+longitude),
         formula=migStatus ~ means + Matern(1|latitude+longitude), 
         HLmethod='ML',data=blackcap)

expect_equal(fl$basicLRT$pvalue,0.008070582,tolerance=1e-5)