### R code from vignette source 'clue.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: clue.Rnw:40-42
###################################################
options(width = 60)
library("clue")


###################################################
### code chunk number 2: clue.Rnw:310-319
###################################################
cl_class_ids.glvq <-
function(x)
    as.cl_class_ids(x$class_ids)
is.cl_partition.glvq <-
function(x)
    TRUE
is.cl_hard_partition.glvq <-
function(x)
    TRUE


###################################################
### code chunk number 3: Cassini-data (eval = FALSE)
###################################################
## data("Cassini")
## plot(Cassini$x, col = as.integer(Cassini$classes),
##      xlab = "", ylab = "")


###################################################
### code chunk number 4: clue.Rnw:889-890
###################################################
data("Cassini")
plot(Cassini$x, col = as.integer(Cassini$classes),
     xlab = "", ylab = "")


###################################################
### code chunk number 5: CKME (eval = FALSE)
###################################################
## data("CKME")
## plot(hclust(cl_dissimilarity(CKME)), labels = FALSE)


###################################################
### code chunk number 6: clue.Rnw:903-904
###################################################
data("CKME")
plot(hclust(cl_dissimilarity(CKME)), labels = FALSE)


###################################################
### code chunk number 7: clue.Rnw:914-916
###################################################
m1 <- cl_medoid(CKME)
table(Medoid = cl_class_ids(m1), "True Classes" = Cassini$classes)


###################################################
### code chunk number 8: Cassini-medoid (eval = FALSE)
###################################################
## plot(Cassini$x, col = cl_class_ids(m1), xlab = "", ylab = "")


###################################################
### code chunk number 9: clue.Rnw:924-925
###################################################
plot(Cassini$x, col = cl_class_ids(m1), xlab = "", ylab = "")


###################################################
### code chunk number 10: clue.Rnw:934-936
###################################################
set.seed(1234)
m2 <- cl_consensus(CKME)


###################################################
### code chunk number 11: clue.Rnw:941-942
###################################################
table(Consensus = cl_class_ids(m2), "True Classes" = Cassini$classes)


###################################################
### code chunk number 12: Cassini-mean (eval = FALSE)
###################################################
## plot(Cassini$x, col = cl_class_ids(m2), xlab = "", ylab = "")


###################################################
### code chunk number 13: clue.Rnw:950-951
###################################################
plot(Cassini$x, col = cl_class_ids(m2), xlab = "", ylab = "")


###################################################
### code chunk number 14: clue.Rnw:984-989
###################################################
data("GVME")
GVME
set.seed(1)
m1 <- cl_consensus(GVME, method = "GV1",
                   control = list(k = 3, verbose = TRUE))


###################################################
### code chunk number 15: clue.Rnw:993-994
###################################################
mean(cl_dissimilarity(GVME, m1, "GV1") ^ 2)


###################################################
### code chunk number 16: clue.Rnw:998-1002
###################################################
data("GVME_Consensus")
m2 <- GVME_Consensus[["MF1/3"]]
mean(cl_dissimilarity(GVME, m2, "GV1") ^ 2)
table(CLUE = cl_class_ids(m1), GV2001 = cl_class_ids(m2))


###################################################
### code chunk number 17: clue.Rnw:1009-1012
###################################################
set.seed(1)
m1 <- cl_consensus(GVME, method = "GV1",
                   control = list(k = 2, verbose = TRUE))


###################################################
### code chunk number 18: clue.Rnw:1016-1019
###################################################
mean(cl_dissimilarity(GVME, m1, "GV1") ^ 2)
m2 <- GVME_Consensus[["MF1/2"]]
mean(cl_dissimilarity(GVME, m2, "GV1") ^ 2)


###################################################
### code chunk number 19: clue.Rnw:1022-1023
###################################################
max(abs(cl_membership(m1) - cl_membership(m2)))


###################################################
### code chunk number 20: clue.Rnw:1027-1029
###################################################
m3 <- cl_consensus(GVME, method = "GV1",
                   control = list(k = 2, verbose = TRUE))


###################################################
### code chunk number 21: clue.Rnw:1032-1033
###################################################
table(GV1 = cl_class_ids(m1), Euclidean = cl_class_ids(m3))


###################################################
### code chunk number 22: clue.Rnw:1036-1037
###################################################
rownames(m1)[cl_class_ids(m1) != cl_class_ids(m3)]


###################################################
### code chunk number 23: clue.Rnw:1061-1066
###################################################
data("Kinship82")
Kinship82
set.seed(1)
m1 <- cl_consensus(Kinship82, method = "GV3",
                   control = list(k = 3, verbose = TRUE))


###################################################
### code chunk number 24: clue.Rnw:1071-1072
###################################################
mean(cl_dissimilarity(Kinship82, m1, "comem") ^ 2)


###################################################
### code chunk number 25: clue.Rnw:1076-1079
###################################################
data("Kinship82_Consensus")
m2 <- Kinship82_Consensus[["JMF"]]
mean(cl_dissimilarity(Kinship82, m2, "comem") ^ 2)


###################################################
### code chunk number 26: clue.Rnw:1083-1085
###################################################
cl_dissimilarity(m1, m2, "comem")
table(CLUE = cl_class_ids(m1), GV2001 = cl_class_ids(m2))


###################################################
### code chunk number 27: clue.Rnw:1088-1089
###################################################
cl_fuzziness(cl_ensemble(m1, m2))


###################################################
### code chunk number 28: clue.Rnw:1109-1111
###################################################
data("Phonemes")
d <- as.dist(1 - Phonemes)


###################################################
### code chunk number 29: clue.Rnw:1115-1116
###################################################
u <- ls_fit_ultrametric(d, control = list(verbose = TRUE))


###################################################
### code chunk number 30: Phonemes (eval = FALSE)
###################################################
## plot(u)


###################################################
### code chunk number 31: clue.Rnw:1126-1127
###################################################
plot(u)


###################################################
### code chunk number 32: clue.Rnw:1137-1138
###################################################
round(cl_dissimilarity(d, u), 4)


###################################################
### code chunk number 33: clue.Rnw:1141-1146
###################################################
hclust_methods <- c("ward", "single", "complete", "average", "mcquitty")
hens <- cl_ensemble(list = lapply(hclust_methods,
                                  function(m) hclust(d, m)))
names(hens) <- hclust_methods
round(sapply(hens, cl_dissimilarity, d), 4)


###################################################
### code chunk number 34: clue.Rnw:1153-1155
###################################################
ahens <- c(L2opt = cl_ensemble(u), hens)
round(cl_dissimilarity(ahens, method = "gamma"), 2)


