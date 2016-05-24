
 library(aster2)

 data(echinacea)

 foo <- echinacea$redata$pop
 bar <- match(levels(foo), as.character(foo))
 baz <- echinacea$redata$id %in% echinacea$redata$id[bar]

 out <- subset(echinacea, baz)

 all(sort(unique(out$redata$id)) == sort(bar))

 opred <- echinacea$repred[baz]
 ipred <- seq(along = echinacea$repred)[baz]
 all(c(0, ipred)[1 + out$repred] == opred)

 opred <- echinacea$regroup[baz]
 all(c(0, ipred)[1 + out$regroup] == opred)

