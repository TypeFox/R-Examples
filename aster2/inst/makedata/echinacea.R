
 ##### this is the script used to make the data set echinacea.rda #####

 library(aster)
 data(echinacea)

 vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
     "hdct02", "hdct03", "hdct04")
 pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
 group <- rep(0, length(pred))
 code <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
 families <- list("bernoulli", "poisson", "zero.truncated.poisson")

 library(aster2, lib.loc = "../../../aster2.Rcheck")

 echinacea <- asterdata(echinacea, vars, pred, group, code, families)

 save(echinacea, file = "echinacea.rda")

