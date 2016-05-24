source("../R/listPartition.r")

cat("\n#######################################################################
########################## Test ListPartition #########################
############################### Creation ###############################
####################################################################\n")

cleanProg(.ListPartition.validity,,,0)
cleanProg(listPartition)
cleanProg(showListPartition)
new("ListPartition")
lcl0 <- lcl3 <- lcl4 <- listPartition()

cleanProg(.listPartition.set,,,2) # CLUSTER_NAMES CRITERION_NAMES
cleanProg(.listPartition.get,,,4) # CLUSTER_NAMES listI length CRITERION_NAMES


lcl3['add'] <- p3a
tryBug(lcl3['sorted'] <- TRUE)
lcl3['add'] <- p3b
lcl3['add'] <- p3c
lcl3['add'] <- p3d
lcl3['add'] <- p3e
lcl3['add'] <- p3f
lcl3['add'] <- p3g
lcl3['add'] <- p3h
lcl3['add'] <- p3i
lcl3['add'] <- p3j

lcl4['add'] <- p4a
lcl4['add'] <- p4b
lcl4['add'] <- p4c
#lcl4['add'] <- p4d
lcl4['add'] <- p4e
lcl4['add'] <- p4f
tryBug(lcl3['c18'] <- c3e)
tryBug(lcl3['ccc18'] <- p3a)


(lcl3)
lcl3['criterionActif'] <- "random"
(lcl3)
tryBug(lcl3['criterionActif'] <- "test")

lcl3['initializationMethod']
lcl3['initializationMethod'] <- "kmeans-"
lcl3['initializationMethod']


lcl3b <- lcl3
lcl3b["c3"] <- "clear"
tryBug(lcl3b["c3"] <- "cleardf")


lcl3['sorted']
tryBug(lcl3['t3'])
lcl3['c9']
tryBug(lcl3[3])
lcl3['random']

lcl3['criterionValues']
lcl3['criterionValuesAsMatrix']
lcl3['criterionValuesAsMatrix',"Calinski.Genolini"]
lcl3['Calinski.Genolini']

lcl3['sorted']


cleanProg(.ListPartition.ordered,,,2) # CRITERION_MIN_OR_MAX length
ordered(lcl3)

lcl3['criterionActif'] <- "Calinski.Harabatz"
ordered(lcl3)
lcl3['criterionActif'] <- "Ray.Turi"
ordered(lcl3)



cleanProg(.ListPartition.plotOne,,,3) # CLUSTER_NAMES length letters
cleanProg(.ListPartition.plotCriterion,,,3) # CLUSTER_NAMES CRITERION_MIN_OR_MAX CRITERION_NAMES
cleanProg(regroup)
lcl3['add'] <- p3e
lcl3['add'] <- p3f
lcl3['add'] <- p3g
lcl3['add'] <- p3h
lcl3['add'] <- p3i
lcl3['add'] <- p3j
lcl3['add'] <- p3e
lcl3['add'] <- p3f
lcl3['add'] <- p3g
lcl3['add'] <- p3h
lcl3['add'] <- p3i
lcl3['add'] <- p3j

(lcl3)
regroup(lcl3)

plot(lcl0)
plot(lcl3)
plot(lcl3,nbCriterion=3,criterion="Calinski.Harabatz")

plot(lcl3)
tryBug(plot(lcl3,criterion=c("random","Calinski.Harabatz")))



par(mfrow=c(2,2))
lcl3['criterionActif'] <- "Calinski.Harabatz"
ordered(lcl3)
plot(lcl3)

lcl3['criterionActif'] <- "random"
ordered(lcl3)
plot(lcl3)

lcl3['criterionActif'] <- "Ray.Turi"
ordered(lcl3)
plot(lcl3)

lcl3['criterionActif'] <- "Davies.Bouldin"
plot(lcl3)

plotCriterion(lcl0)
plotCriterion(lcl3)
ordered(lcl3)
plotCriterion(lcl3,criterion=c("Calinski.Harabatz","Calinski.Genolini","Davies.Bouldin"))



cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++ Fin Test ListPartition +++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

