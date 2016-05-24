################################################################################
##
## $Id: map.market.test.R 346 2006-10-01 05:08:55Z enos $
##
## Test for the map.market function
##
################################################################################

library(portfolio)

load("map.market.test.RData")

## save(test.data, truth, file = "map.market.test.RData", compress = TRUE)

result <- map.market(id = test.data$id, area = test.data$area,
                     group = test.data$group, color = test.data$color,
                     lab = TRUE, print = FALSE)

stopifnot(
          isTRUE(all.equal(getGrob(result, "MAP"), truth))
          )
