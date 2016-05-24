### R code from vignette source 'QuACN.Rnw'

###################################################
### code chunk number 1: QuACN.Rnw:64-65 (eval = FALSE)
###################################################
## install.packages("QuACN")


###################################################
### code chunk number 2: QuACN.Rnw:72-74
###################################################
library("graph")
library("RBGL")


###################################################
### code chunk number 3: QuACN.Rnw:78-79
###################################################
library("QuACN")


###################################################
### code chunk number 4: randGraph
###################################################

set.seed(666)
g <- randomGraph(1:8, 1:5, 0.36, weights=FALSE)
g


###################################################
### code chunk number 5: QuACN.Rnw:90-91
###################################################
A <- adjacencyMatrix(g)


###################################################
### code chunk number 6: randGraph
###################################################
A
g <- as(A, "graphNEL")
g


###################################################
### code chunk number 7: QuACN.Rnw:111-115
###################################################
nodeDataDefaults(g, "atom") <- "C"
nodeData(g, "6", "atom") <- "O"
edgeDataDefaults(g, "bond") <- 1
edgeData(g, "2", "3", "bond") <- 2


###################################################
### code chunk number 8: QuACN.Rnw:130-133
###################################################
g2 <- randomGraph(paste("A", 1:100, sep=""), 1:4, p=0.03, weights=FALSE)
lcc <- getLargestSubgraph(g2)
lcc


###################################################
### code chunk number 9: QuACN.Rnw:143-145
###################################################
sg.1ed <- edgeDeletedSubgraphs(g)
sg.2ed <- edgeDeletedSubgraphs(sg.1ed)


###################################################
### code chunk number 10: QuACN.Rnw:158-163
###################################################
mat.adj <- adjacencyMatrix(g)
mat.dist <- distanceMatrix(g)
vec.degree <- graph::degree(g)
ska.dia <- diameter(g)
ska.dia <- diameter(g, mat.dist)


###################################################
### code chunk number 11: QuACN.Rnw:182-184
###################################################
wien <- wiener(g)
wiener(g, mat.dist)


###################################################
### code chunk number 12: QuACN.Rnw:192-194
###################################################
harary(g)
harary(g, mat.dist)


###################################################
### code chunk number 13: QuACN.Rnw:201-203
###################################################
balabanJ(g)
balabanJ(g, mat.dist)


###################################################
### code chunk number 14: QuACN.Rnw:234-237
###################################################
compactness(g)
compactness(g, mat.dist)
compactness(g, mat.dist, wiener(g, mat.dist))


###################################################
### code chunk number 15: QuACN.Rnw:246-250
###################################################
productOfRowSums(g, log=FALSE)
productOfRowSums(g, log=TRUE)
productOfRowSums(g, mat.dist, log=FALSE)
productOfRowSums(g, mat.dist, log=TRUE)


###################################################
### code chunk number 16: QuACN.Rnw:258-261
###################################################
hyperDistancePathIndex(g)
hyperDistancePathIndex(g, mat.dist)
hyperDistancePathIndex(g, mat.dist, wiener(g, mat.dist))


###################################################
### code chunk number 17: QuACN.Rnw:273-276
###################################################
dob <- dobrynin(g)
dob <- dobrynin(g, mat.dist)
dob$eccentricityVertex


###################################################
### code chunk number 18: QuACN.Rnw:286-287
###################################################
dob$eccentricityGraph


###################################################
### code chunk number 19: QuACN.Rnw:297-298
###################################################
dob$avgeccOfG


###################################################
### code chunk number 20: QuACN.Rnw:307-308
###################################################
dob$ecentricVertex


###################################################
### code chunk number 21: QuACN.Rnw:318-319
###################################################
dob$ecentricGraph


###################################################
### code chunk number 22: QuACN.Rnw:328-329
###################################################
dob$vertexCentrality


###################################################
### code chunk number 23: QuACN.Rnw:339-340
###################################################
dob$graphIntegration


###################################################
### code chunk number 24: QuACN.Rnw:350-351
###################################################
dob$unipolarity


###################################################
### code chunk number 25: QuACN.Rnw:360-361
###################################################
dob$vertexDeviation


###################################################
### code chunk number 26: QuACN.Rnw:371-372
###################################################
dob$variation


###################################################
### code chunk number 27: QuACN.Rnw:382-383
###################################################
dob$centralization


###################################################
### code chunk number 28: QuACN.Rnw:393-394
###################################################
dob$avgDistance


###################################################
### code chunk number 29: QuACN.Rnw:403-404
###################################################
dob$distVertexDeviation


###################################################
### code chunk number 30: QuACN.Rnw:414-415
###################################################
dob$meanDistVertexDeviation


###################################################
### code chunk number 31: QuACN.Rnw:432-434
###################################################
totalAdjacency(g)
totalAdjacency(g, mat.adj)


###################################################
### code chunk number 32: QuACN.Rnw:456-466
###################################################
zagreb1(g)
zagreb1(g, vec.degree)
zagreb2(g)
zagreb2(g, vec.degree)
modifiedZagreb(g)
modifiedZagreb(g, vec.degree)
augmentedZagreb(g)
augmentedZagreb(g, vec.degree)
variableZagreb(g)
variableZagreb(g, vec.degree)


###################################################
### code chunk number 33: QuACN.Rnw:474-476
###################################################
randic(g)
randic(g, vec.degree)


###################################################
### code chunk number 34: QuACN.Rnw:486-489
###################################################
complexityIndexB(g)
complexityIndexB(g, mat.dist)
complexityIndexB(g, mat.dist, vec.degree)


###################################################
### code chunk number 35: QuACN.Rnw:496-498
###################################################
normalizedEdgeComplexity(g)
normalizedEdgeComplexity(g, totalAdjacency(g, mat.adj))


###################################################
### code chunk number 36: QuACN.Rnw:506-508
###################################################
atomBondConnectivity(g)
atomBondConnectivity(g, vec.degree)


###################################################
### code chunk number 37: QuACN.Rnw:529-535
###################################################
geometricArithmetic1(g)
geometricArithmetic1(g, vec.degree)
geometricArithmetic2(g)
geometricArithmetic2(g, mat.dist)
geometricArithmetic3(g)
geometricArithmetic3(g, mat.dist)


###################################################
### code chunk number 38: QuACN.Rnw:543-545
###################################################
narumiKatayama(g)
narumiKatayama(g, vec.degree)


###################################################
### code chunk number 39: QuACN.Rnw:561-564
###################################################
topologicalInfoContent(g)
topologicalInfoContent(g, mat.dist)
topologicalInfoContent(g, mat.dist, vec.degree)


###################################################
### code chunk number 40: QuACN.Rnw:578-588
###################################################
#I_D(G)
bonchev1(g)
bonchev1(g, mat.dist)
#I^W_D(G)
bonchev2(g)
bonchev2(g, mat.dist)
bonchev2(g, mat.dist, wiener(g))
#I^E_D(G)
bonchev3(g)
bonchev3(g, mat.dist)


###################################################
### code chunk number 41: QuACN.Rnw:598-601
###################################################
bertz(g)
bertz(g, mat.dist)
bertz(g, mat.dist, vec.degree)


###################################################
### code chunk number 42: QuACN.Rnw:610-612
###################################################
radialCentric(g)
radialCentric(g, mat.dist)


###################################################
### code chunk number 43: QuACN.Rnw:621-623
###################################################
vertexDegree(g)
vertexDegree(g, vec.degree)


###################################################
### code chunk number 44: QuACN.Rnw:645-651
###################################################
#Balaban-like information index U(G)
balabanlike1(g)
balabanlike1(g, mat.dist)
#Balaban-like information index X(G)
balabanlike2(g)
balabanlike2(g, mat.dist)


###################################################
### code chunk number 45: QuACN.Rnw:664-666
###################################################
graphVertexComplexity(g)
graphVertexComplexity(g, mat.dist)


###################################################
### code chunk number 46: QuACN.Rnw:680-682
###################################################
graphDistanceComplexity(g)
graphDistanceComplexity(g, mat.dist)


###################################################
### code chunk number 47: QuACN.Rnw:694-695
###################################################
informationBondIndex(g)


###################################################
### code chunk number 48: QuACN.Rnw:705-707
###################################################
edgeEqualityMIC(g)
edgeEqualityMIC(g, vec.degree)


###################################################
### code chunk number 49: QuACN.Rnw:719-721
###################################################
edgeMagnitudeMIC(g)
edgeMagnitudeMIC(g, vec.degree)


###################################################
### code chunk number 50: QuACN.Rnw:731-732
###################################################
symmetryIndex(g)


###################################################
### code chunk number 51: QuACN.Rnw:742-744
###################################################
distanceDegreeMIC(g)
distanceDegreeMIC(g, mat.dist)


###################################################
### code chunk number 52: QuACN.Rnw:754-756
###################################################
distanceDegreeEquality(g)
distanceDegreeEquality(g, mat.dist)


###################################################
### code chunk number 53: QuACN.Rnw:766-768
###################################################
distanceDegreeCompactness(g)
distanceDegreeCompactness(g, mat.dist)


###################################################
### code chunk number 54: QuACN.Rnw:778-780
###################################################
informationLayerIndex(g)
informationLayerIndex(g, mat.dist)


###################################################
### code chunk number 55: QuACN.Rnw:811-812
###################################################
mediumArticulation(g)


###################################################
### code chunk number 56: QuACN.Rnw:824-826
###################################################
efficiency(g)
efficiency(g, mat.dist)


###################################################
### code chunk number 57: QuACN.Rnw:837-838
###################################################
graphIndexComplexity(g)


###################################################
### code chunk number 58: QuACN.Rnw:855-857
###################################################
offdiagonal(g)
offdiagonal(g, vec.degree)


###################################################
### code chunk number 59: QuACN.Rnw:882-884
###################################################
spanningTreeSensitivity(g)
spanningTreeSensitivity(g, sg.1ed)


###################################################
### code chunk number 60: QuACN.Rnw:900-904
###################################################
distanceDegreeCentric(g)
distanceDegreeCentric(g, mat.dist)
distanceCodeCentric(g)
distanceCodeCentric(g, mat.dist)


###################################################
### code chunk number 61: QuACN.Rnw:1014-1021
###################################################
l1 <- infoTheoreticGCM(g)
l2 <- infoTheoreticGCM(g, mat.dist, coeff="lin", infofunct="sphere", lambda=1000)
l3 <- infoTheoreticGCM(g, mat.dist, coeff="const", infofunct="pathlength", lambda=4000)
l4 <- infoTheoreticGCM(g, mat.dist, coeff="quad", infofunct="vertcent", lambda=1000)
l5 <- infoTheoreticGCM(g, mat.dist, coeff="exp", infofunct="degree", lambda=1000)
l1
l5


###################################################
### code chunk number 62: QuACN.Rnw:1026-1030
###################################################
l5mpfr <- infoTheoreticGCM(g, mat.dist, coeff="exp", infofunct="degree", lambda=1000, prec=128)
l5mpfr$entropy
l5mpfr$entropy * 2^3
as.double(l5mpfr$entropy * 2^3)


###################################################
### code chunk number 63: QuACN.Rnw:1067-1073
###################################################
lv1 <- infoTheoreticLabeledV1(g, coeff="exp")
lv1$entropy
lv2 <- infoTheoreticLabeledV2(g, ci=list(`C` = 0.8, `O` = 1))
lv2$entropy
le <- infoTheoreticLabeledE(g, coeff="quad")
le$entropy


###################################################
### code chunk number 64: QuACN.Rnw:1078-1082
###################################################
lv1e <- infoTheoreticSum(lv1, le)
lv1e$entropy
lv2e <- infoTheoreticSum(lv2, le)
lv2e$entropy


###################################################
### code chunk number 65: QuACN.Rnw:1114-1115
###################################################
eigenvalueBased(g, adjacencyMatrix,2)


###################################################
### code chunk number 66: QuACN.Rnw:1118-1119
###################################################
eigenvalueBased(g, laplaceMatrix,2)


###################################################
### code chunk number 67: QuACN.Rnw:1122-1123
###################################################
eigenvalueBased(g, distanceMatrix,2)


###################################################
### code chunk number 68: QuACN.Rnw:1126-1127
###################################################
eigenvalueBased(g,distancePathMatrix,2)


###################################################
### code chunk number 69: QuACN.Rnw:1130-1131
###################################################
eigenvalueBased(g, augmentedMatrix,2)


###################################################
### code chunk number 70: QuACN.Rnw:1134-1135
###################################################
eigenvalueBased(g, extendedAdjacencyMatrix,2)


###################################################
### code chunk number 71: QuACN.Rnw:1138-1139
###################################################
eigenvalueBased(g, vertConnectMatrix,2) 


###################################################
### code chunk number 72: QuACN.Rnw:1142-1143
###################################################
eigenvalueBased(g, randomWalkMatrix,2)  


###################################################
### code chunk number 73: QuACN.Rnw:1146-1147
###################################################
eigenvalueBased(g, weightStrucFuncMatrix_lin,2) 


###################################################
### code chunk number 74: QuACN.Rnw:1150-1151
###################################################
eigenvalueBased(g, weightStrucFuncMatrix_exp,2)


###################################################
### code chunk number 75: QuACN.Rnw:1167-1169
###################################################
energy(g)
laplacianEnergy(g)


###################################################
### code chunk number 76: QuACN.Rnw:1180-1182
###################################################
estrada(g)
laplacianEstrada(g)


###################################################
### code chunk number 77: QuACN.Rnw:1190-1191
###################################################
spectralRadius(g)


###################################################
### code chunk number 78: QuACN.Rnw:1210-1212
###################################################
oneEdgeDeletedSubgraphComplexity(g)
oneEdgeDeletedSubgraphComplexity(g, sg.1ed)


###################################################
### code chunk number 79: QuACN.Rnw:1225-1227
###################################################
twoEdgesDeletedSubgraphComplexity(g)
twoEdgesDeletedSubgraphComplexity(g, sg.2ed)


###################################################
### code chunk number 80: QuACN.Rnw:1236-1238
###################################################
localClusteringCoeff(g)
localClusteringCoeff(g, deg=vec.degree)


###################################################
### code chunk number 81: QuACN.Rnw:1247-1250
###################################################
loccc <- localClusteringCoeff(g)
globalClusteringCoeff(g)
globalClusteringCoeff(g, loc=loccc)


###################################################
### code chunk number 82: QuACN.Rnw:1270-1272
###################################################
connectivityID(g)
connectivityID(g, deg=vec.degree)


###################################################
### code chunk number 83: QuACN.Rnw:1283-1285
###################################################
minConnectivityID(g)
minConnectivityID(g, deg=vec.degree)


###################################################
### code chunk number 84: QuACN.Rnw:1300-1302
###################################################
primeID(g)
primeID(g, deg=vec.degree)


###################################################
### code chunk number 85: QuACN.Rnw:1319-1320
###################################################
bondOrderID(g)


###################################################
### code chunk number 86: QuACN.Rnw:1336-1338
###################################################
balabanID(g)
balabanID(g, dist=mat.dist)


###################################################
### code chunk number 87: QuACN.Rnw:1349-1351
###################################################
minBalabanID(g)
minBalabanID(g, dist=mat.dist)


###################################################
### code chunk number 88: QuACN.Rnw:1367-1368
###################################################
weightedID(g)


###################################################
### code chunk number 89: QuACN.Rnw:1388-1390
###################################################
huXuID(g)
huXuID(g, deg=vec.degree)


###################################################
### code chunk number 90: QuACN.Rnw:1408-1411
###################################################
calculateDescriptors(g, "wiener")
calculateDescriptors(g, 1001)
calculateDescriptors(g, 2000, labels=TRUE)


###################################################
### code chunk number 91: QuACN.Rnw:1415-1416
###################################################
sessionInfo()


