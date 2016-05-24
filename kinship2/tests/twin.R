require(kinship2)
#
# Test some twins data from Curtis Oswold
#
twindat <- c(1,3,4,2,
             2,0,0,1,
             3,8,7,1,
             4,6,5,2,
             5,0,0,2,
             6,0,0,1,
             7,0,0,2,
             8,0,0,1,
             100,3,4,1,
             101,3,4,2,
             102,3,4,2,
             103,3,4,2,
             104,3,4,2,
             105,3,4,2,
             106,3,4,2,
             107,0,0,1,
             108,0,0,1,
             201,2,1,1,
             202,2,1,1,
             203,2,1,1,
             204,2,1,1,
             205,107,102,1,
             206,108,103,2)
twindat <- matrix(twindat, ncol=4, byrow=T)
dimnames(twindat) <- list(NULL, c('id', 'dadid', 'momid', 'sex'))
twindat <- data.frame(twindat)

## set up a fraternal twin set, and a set of triplets with kids from
## their marriages to test kinship coeff
relate=data.frame(id1=c(101,102,104,203), id2=c(102,103,105,204), code=c(1,1,2,1))

tped <- with(twindat, pedigree(id, dadid, momid, sex,
                               relation=relate))


## plot(tped)

## should show kinship coeff of 0.5 for where MZ twins are
## ids: 102-103 and 203-204
kinmat <- kinship(tped)

kinmat[c(10:16,19:23),c(10:16,19:23)]
                                       



## simple test case for kinship of MZ twins from Claus Ekstrom, 9/2012
mydata <- data.frame(id=1:4, dadid=c(NA, NA, 1, 1),
momid=c(NA, NA, 2, 2), sex=c("male", "female", "male", "male"),
famid=c(1,1,1,1))
relation <- data.frame(id1=c(3), id2=c(4), famid=c(1), code=c(1))

x <- pedigree(id=mydata$id, dadid=mydata$dadid, momid=mydata$momid, sex=mydata$sex, relation=relation)

#plot(x)

kinout <- kinship(x)
kinship2:::kinship.pedigree(x)

