#######################
## Tests on PlabStat ##
#######################


load("unrepTrit.RData")

load("unrepTritPlabStat.RData")

unrepTrit <- data.matrix(unrepTrit)

fieldMap <- matrix(nrow = max(unrepTrit[,1]),
                   ncol = max(unrepTrit[,2]))

fieldMap[unrepTrit[,1:2]] <- unrepTrit[,4]


shape <- list(c(1),
              c(1),
              c(1:4),
              c(1:4))

library(mvngGrAd)

testUnrepTrit <- movingGrid(rows=unrepTrit[,1],
                            columns = unrepTrit[,2],
                            obsPhe = unrepTrit[,4],
                            shapeCross = shape,
                            layers = NULL)



all.equal(fitted(testUnrepTrit),unrepTritPlabStat[,4],
          check.attributes = FALSE,
          check.names = FALSE)





