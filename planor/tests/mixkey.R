## Asymmetrical (mixed) factorial design
## cf Planor Manual, Example 2 page 12
library("planor")
# Four treatment factors at 6, 6, 4, 2 levels and one 6-level block factor
# Model: block+(A+B+C+D)^2 ; Estimate: A+B+C+D\n")
# N = 144 = 2^4 x 3^2 experimental units

mixKey <- planor.designkey(factors=c( LETTERS[1:4], "block"), 
                           nlevels=c(6,6,4,2,6), 
                           block=~block,
                           model=~block+(A+B+C+D)^2, 
                           estimate=~A+B+C+D, 
                           nunits=144,
                       base=~A+B+D, max.sol=2)

## Tests on the listofkeyrings class
summary(mixKey)
alias(mixKey)
mixPlan <- planor.design(key=mixKey, select=c(1,1), randomize=~block/UNITS)
# resultats aleatoires: print(getDesign(mixPlan)[1:25,])

## Tests on the designkey class
summary(mixKey[c(1, 1)])
alias(mixKey[c(1, 1)])
mixPlan <- planor.design(key=mixKey, select=c(1,1), randomize=~block/UNITS)
# resultats aleatoires: print(getDesign(mixPlan)[1:25,])
