library(RXMCDA)

performanceTable <- rbind(c(1,2,3),c(4,5,6))

rownames(performanceTable) <- c("x","y")

colnames(performanceTable) <- c("g1","g2","g3")

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putPerformanceTable(tree,performanceTable)

perfTable2 <- getPerformanceTables(tree)

stopifnot(identical(performanceTable,perfTable2[[1]]))
