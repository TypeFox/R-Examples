library("planor")
#---------------------------------------------------------------------------
# Aucun terme ineligible
#---------------------------------------------------------------------------
F0 <- planor.factors( factors=c(LETTERS[1:5], "bloc"), nlevels=rep(2,6) )
K0 <- planor.designkey(factors=c(LETTERS[1:5], "bloc"), nlevels=rep(2,6),
                       model=~bloc+(A+B+C+D+E)^2, estimate=~A+B+C ,
                       nunits=2^3,
                       base=~A+B+C, max.sol=2, verbose=T)
