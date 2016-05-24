### R code from vignette source 'symmoments.rnw'

###################################################
### code chunk number 1: lib
###################################################
library(symmoments)


###################################################
### code chunk number 2: m1234
###################################################
writeLines(toLatex(callmultmoments(c(1,2,3,4))))


###################################################
### code chunk number 3: callmultmoments
###################################################
m1234 <- callmultmoments(c(1,2,3,4))
unclass(m1234)


###################################################
### code chunk number 4: print-method
###################################################
m1234


###################################################
### code chunk number 5: toLatex-method
###################################################
toLatex(m1234)


###################################################
### code chunk number 6: toLatex-method2 (eval = FALSE)
###################################################
## writeLines(toLatex(m1234), "yourfilename")


###################################################
### code chunk number 7: covariance
###################################################
matrix(c(4,2,1,1,2,3,1,1,1,1,2,1,1,1,1,2),nrow=4)


###################################################
### code chunk number 8: evaluate-method
###################################################
evaluate(m1234, c(4, 2, 1, 1, 3, 1, 1, 2, 1, 2))


###################################################
### code chunk number 9: simulate-method
###################################################
#simulate(m1234,1000,NULL, c(1,2,0,3),c(4,2,1,1,2,3,1,1,1,1,2,1,1,1,1,2))


###################################################
### code chunk number 10: noncentral1
###################################################
# as.matrix(toLatex_noncentral(c(1,3))) 


###################################################
### code chunk number 11: noncentral2
###################################################
# make.all.moments(c(3,3))


###################################################
### code chunk number 12: noncentral2a
###################################################
# evaluate_noncentral(c(1,3),c(1,2),c(1,0,1))


###################################################
### code chunk number 13: noncentral3a
###################################################
library(mpoly)  # required for this example


###################################################
### code chunk number 14: noncentral3
###################################################
t0 <- mpoly(list(c(coef=3,x1=2),c(coef=2,x1=1,x2=3),
                    c(coef=-4,z=2),c(coef=1,x1=1,x2=2,z=1)))                       
print(t0)


###################################################
### code chunk number 15: <noncentral4
###################################################
t1 <<- convert.mpoly(t0)      
t1


###################################################
### code chunk number 16: noncentral5
###################################################
t2 <<- convert.multipol(t1)      
t2


###################################################
### code chunk number 17: noncentral6
###################################################
print(mpoly(convert.mpoly(convert.multipol(t2))))     


###################################################
### code chunk number 18: noncentral7
###################################################
# evaluate_expected.polynomial(t0,c(1,2,3),c(1,0,0,1,0,1))


###################################################
### code chunk number 19: trees1
###################################################
exam.Newick      <- "(((a,b),c),d);"


###################################################
### code chunk number 20: trees2
###################################################
library(ape)
exam.phylo       <- read.tree(text=exam.Newick)
exam.phylo


###################################################
### code chunk number 21: trees3
###################################################
exam.matching    <- as.matching(exam.phylo)
exam.matching 


###################################################
### code chunk number 22: trees3
###################################################
exam.L.matrix    <- toMoment(exam.matching)
exam.L.matrix


###################################################
### code chunk number 23: trees3
###################################################
# backto.matching  <- toMatching(exam.L.matrix )
# backto.L.matrix  <- toMoment(backto.matching$matching[,c(1,2)])
# backto.Newick  <- toNewick(exam.L.matrix$L,type="square")
# backto.L.matrix.tips <- toMoment(exam.Newick,tip.label=c("d","a","c","b"))
# backto.phylo     <- as.phylo.matching(backto.matching)


