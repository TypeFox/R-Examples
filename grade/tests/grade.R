library(grade)

set.seed(42)

x <- c(1,2)
sans <- "[1,2]"
grade.interval(x, sans)

grade.interval(c(-Inf, Inf), "[-Inf, Inf]", useinf=T)

grade.negative(NULL, -Inf)
grade.negative(NULL, 0)
grade.negative(NULL, -1)
grade.negative(NULL, -Inf, useinf=T)

grade.number(NA, "NA", usena=T)
grade.number(1, 1.1, tolerance=.01)
grade.number(1, 1.1, tolerance=.100001)

x1 <- runif(1)
x <- c(x1, 1-x1)
x
grade.discreteprobability(NULL, x, checkcorrect=F)
grade.discreteprobability(x, c(x[2], x[1]), ordered=T)
grade.discreteprobability(x, c(x[2], x[1]), ordered=F)

grade.set(c(1,2,3,4), c(2,4,3,1))

grade.orderedset(c(1,2,3,4), "[1,2,3,4]")
grade.orderedset(c(1,2,3,4), "[1,3,2,4]")
