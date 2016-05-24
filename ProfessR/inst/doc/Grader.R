### R code from vignette source 'Grader.Rnw'

###################################################
### code chunk number 1: Grader.Rnw:83-86
###################################################
g = rnorm(n=200, m=82, sd=10)
g[g>100] = 100
g[g<1] = 1


###################################################
### code chunk number 2: Grader.Rnw:94-96
###################################################
B = boxplot(g, plot=FALSE)
divs = c(min(g), B$stats[1:4] + diff(B$stats)/2, max(g) )


###################################################
### code chunk number 3: Grader.Rnw:104-111
###################################################

library(ProfessR)
## get(getOption("device"))(width = 12, height = 7)
dev.new()


D1 = do.grades(g, divs=divs, tit="GEOL 105 Exam 1")


###################################################
### code chunk number 4: Grader.Rnw:155-170
###################################################

data(E2grades)

g = E2grades

B = boxplot(g[g>1], plot=FALSE)
divs = c(min(g), B$stats[1:4] + diff(B$stats)/2, max(g) )
##   get(getOption("device"))(width = 12, height = 7)


G1 = do.grades(g, divs=divs , tit="GEOL 105 Exam 1")


J = jist(G1$hist, G1$grades, G1$lett)



