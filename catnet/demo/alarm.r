library(catnet)

data(alarmnet)
cnPlot(alarmnet)

alarm.data <- cnSamples(alarmnet, 10000)
cnLoglik(alarmnet, alarm.data)

eval <- cnSearchOrder(data=alarm.data,perturbations=NULL, maxParentSet=4, parentSizes=NULL, maxComplexity=2*cnComplexity(alarmnet), nodeOrder=cnOrder(alarmnet), echo=TRUE)

bstnet<- cnFind(eval, cnComplexity(alarmnet))
bstnet
cnCompare(alarmnet, bstnet)

cnCompare(alarmnet, cnFindBIC(eval))

ecmp <- cnCompare(alarmnet, eval,extended=TRUE)
pdf("alarmeval.pdf")
cnPlot(ecmp)
dev.off()

cnSetSeed(123)
sanets1 <- cnSearchSA(data = alarm.data, perturbations = NULL, maxParentSet = 2, 
     maxComplexity = 600, parentsPool = NULL, fixedParents = NULL, 
     selectMode = "BIC", tempStart = 0, tempCoolFact = 1, tempCheckOrders = 1000, 
     maxIter = 1000, orderShuffles = 0, stopDiff = 1e-10, numThreads=4, priorSearch = NULL, echo=TRUE)
sanets1

sanet1 <- cnFindBIC(sanets1)
cnCompare(alarmnet, sanet1)

sanet1@meta <- "sanet1, SA stage I"
sanet1
cnPlot(sanet1, "sanet1")

cnSetSeed(123)
sanets <- cnSearchSA(data = alarm.data, perturbations = NULL, maxParentSet = 3, 
     maxComplexity = 600, parentsPool = NULL, fixedParents = NULL, 
     selectMode = "BIC", tempStart = 1e-2, tempCoolFact = 0.9, 
     tempCheckOrders = 20, maxIter = 200, orderShuffles = 4, stopDiff = 1e-2, 
     numThreads=4, priorSearch = sanets1, echo=TRUE)
sanets

sanet2 <- cnFindBIC(sanets)
cnCompare(alarmnet, sanet2)

sanet2@meta <- "sanet2, SA stage II"
sanet2
cnPlot(sanet2, "sanet2")

