context("test primarySuppression()")

sp <- searchpaths()
fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
problem <- get(load(file=fn))

p1 <- primarySuppression(problem, type="freq", maxN=2)
p1.sdc <- get.problemInstance(p1@problemInstance, "sdcStatus")
expect_that(p1, is_a("sdcProblem"))
expect_that(sum(p1.sdc=="u"), equals(1))
expect_that(which(p1.sdc=="u"), equals(6))

expect_that(primarySuppression(problem, type="nk", n=2, k=80), throws_error())
p2 <- primarySuppression(problem, type="nk", n=2, k=75, numVarInd=1)
p2.sdc <- get.problemInstance(p2@problemInstance, "sdcStatus")
expect_that(p2, is_a("sdcProblem"))
expect_that(sum(p2.sdc=="u"), equals(2))
expect_that(which(p2.sdc=="u"), equals(c(6,14)))

expect_that(primarySuppression(problem, type="p", p=80), throws_error())
expect_that(primarySuppression(problem, type="p", p=0, numVarInd=1), throws_error())
expect_that(primarySuppression(problem, type="p", p=100, numVarInd=1), throws_error())
p3 <- primarySuppression(problem, type="p", p=70, numVarInd=1)
p3.sdc <- get.problemInstance(p3@problemInstance, "sdcStatus")
expect_that(p3, is_a("sdcProblem"))
expect_that(sum(p3.sdc=="u"), equals(2))
expect_that(which(p3.sdc=="u"), equals(c(6,14)))

expect_that(primarySuppression(problem, type="pq", pq=c(80,90)), throws_error())
expect_that(primarySuppression(problem, type="pq", pq=c(85, 80), numVarInd=1), throws_error())
expect_that(primarySuppression(problem, type="pq", pq=c(85,180), numVarInd=1), throws_error())
expect_that(primarySuppression(problem, type="pq", pq=c(110, 120), numVarInd=1), throws_error())
p4 <- primarySuppression(problem, type="pq", pq=c(60, 80), numVarInd=1)
p4.sdc <- get.problemInstance(p4@problemInstance, "sdcStatus")
expect_that(p4, is_a("sdcProblem"))
expect_that(sum(p4.sdc=="u"), equals(2))
expect_that(which(p4.sdc=="u"), equals(c(6,14)))

problem@dataObj@rawData$val[5] <- -5
expect_that(primarySuppression(problem, type="nk", n=2, k=80, numVarInd=1), throws_error())
expect_that(primarySuppression(problem, type="p", p=70, numVarInd=1), throws_error())
expect_that(primarySuppression(problem, type="pq", pq=c(60,80), numVarInd=1), throws_error())


