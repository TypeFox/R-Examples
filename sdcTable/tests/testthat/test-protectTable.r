context("test protectTable()")

sp <- searchpaths()
fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
problem <- get(load(file=fn))

p1 <- protectTable(problem, method="OPT", useC=FALSE, verbose=FALSE)
p2 <- protectTable(problem, method="OPT", useC=TRUE, verbose=FALSE)
expect_that(p1@finalData, is_equivalent_to(p2@finalData))
expect_that(p1, is_a("safeObj"))
expect_that(p1@nrPublishableCells, equals(11))
expect_that(p1@nrSecondSupps, equals(3))
expect_that(which(p1@finalData$sdcStatus!="s"), equals(c(5,6,11,12)))

expect_that(p2, is_a("safeObj"))
expect_that(p2@nrPublishableCells, equals(11))
expect_that(p2@nrSecondSupps, equals(3))
expect_that(which(p2@finalData$sdcStatus!="s"), equals(c(5,6,11,12)))

p3 <- protectTable(problem, method="HYPERCUBE", verbose=FALSE)
expect_that(which(p3@finalData$sdcStatus!="s"), equals(c(5,6,11,12)))

p4 <- protectTable(problem, method="HITAS", useC=TRUE, verbose=FALSE)
expect_that(which(p4@finalData$sdcStatus!="s"), equals(c(5,6,11,12)))

rm(problem)
fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
problem <- get(load(file=fn))
p5 <- protectTable(problem, method="OPT", useC=TRUE, verbose=FALSE)
expect_that(p5, is_a("safeObj"))
expect_that(p5@nrPublishableCells, equals(15))
expect_that(p5@nrSecondSupps, equals(0))
expect_that(sum(p5@finalData$sdcStatus!="s"), equals(0))






