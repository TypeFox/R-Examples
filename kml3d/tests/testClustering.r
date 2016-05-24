### clusterization est une partition associé a une longData, ou une clusterizLongData.
### cet objet ne devrait pouvoir exister que dans un cld

source("testClustering.data.r")

cat("####################################################################
########################## Test  Clustering ########################
############################## Creation ############################
####################################################################\n")
cleanProg(calculCriterion)
calculCriterion(ld2,p2a)
calculCriterion(ld2,p2a,criterionNames="calinski")
calculCriterion(ld2,p2a,criterionNames="test")

calculCriterion(ld3n,resizePartition(ld3n,p3a))
calculCriterion(ld3n,resizePartition(ld3n,p3a),imputationMethod="LOCF")
calculCriterion(ld3n,resizePartition(ld3n,p3a),imputationMethod="LOCB")
calculCriterion(ld3n,resizePartition(ld3n,p3a),imputationMethod="LI-LOCBF")


cleanProg(.Clustering.validity,,,)
new("Clustering")



cat("\n####################################################################
########################## Test  Clustering ########################
############################# Accesseurs ###########################
####################################################################\n")

# Héritage
cleanProg(.clustering.get)
cleanProg(.clustering.set)
c2an["nbClusters"]
c2a["clusters"]
c3d["nbClusters"]
c3en["clusters"]
c2an["calinski"]
c2a["percentEachCluster"]
c2an["convergenceTime"]
c2an["multiplicity"]
c2an["criterionNames"]
c2an["criterionValues"]
c2an["algorithm"]

try(c2an["calinksi"])
c3an["calinski"]
c3an["test"]
c3an["percentEachCluster"]
c3an["convergenceTime"]
c3an["criterionNames"]
c3an["criterionValues"]
c3an["algorithm"]

c3an["convergenceTime"] <- 5
c3an["multiplicity"] <- c(c3an["multiplicity"]+2)
try(c3an["criterionNames"] <- c(c3an["criterionNames"],"test2"))
try(c3an["criterionValues"] <- c(c3an["criterionValues"],test2=rnorm(1,2)))

cat("####################################################################
########################### Test Clustering ########################
############################## Affichage ###########################
####################################################################\n")

cleanProg(.Clustering.show,,,1) #LETTERS
c2an
c4d

cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++++++ Fin Test Clustering ++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

