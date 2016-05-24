### R code from vignette source 'recommenderlab.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: recommenderlab.Rnw:101-104
###################################################
options(scipen=3, digits=4, prompt="R> ", eps=FALSE, width=75)
### for sampling
set.seed(1234)


###################################################
### code chunk number 2: recommenderlab.Rnw:1058-1059
###################################################
library("recommenderlab")


###################################################
### code chunk number 3: recommenderlab.Rnw:1067-1072
###################################################
m <- matrix(sample(c(as.numeric(0:5), NA), 50, 
    replace=TRUE, prob=c(rep(.4/6,6),.6)), ncol=10,
    dimnames=list(user=paste("u", 1:5, sep=''), 
	item=paste("i", 1:10, sep='')))
m


###################################################
### code chunk number 4: recommenderlab.Rnw:1079-1082
###################################################
r <- as(m, "realRatingMatrix")
r
#as(r,"dgCMatrix")


###################################################
### code chunk number 5: recommenderlab.Rnw:1087-1088
###################################################
identical(as(r, "matrix"),m)


###################################################
### code chunk number 6: recommenderlab.Rnw:1096-1098
###################################################
as(r, "list")
head(as(r, "data.frame"))


###################################################
### code chunk number 7: recommenderlab.Rnw:1111-1113
###################################################
r_m <- normalize(r)
r_m


###################################################
### code chunk number 8: recommenderlab.Rnw:1119-1121 (eval = FALSE)
###################################################
## image(r, main = "Raw Ratings")
## image(r_m, main = "Normalized Ratings")


###################################################
### code chunk number 9: image1
###################################################
print(image(r, main = "Raw Ratings"))


###################################################
### code chunk number 10: image2
###################################################
print(image(r_m, main = "Normalized Ratings"))


###################################################
### code chunk number 11: recommenderlab.Rnw:1154-1156
###################################################
r_b <- binarize(r, minRating=4)
as(r_b, "matrix")


###################################################
### code chunk number 12: recommenderlab.Rnw:1168-1170
###################################################
data(Jester5k)
Jester5k


###################################################
### code chunk number 13: recommenderlab.Rnw:1178-1180
###################################################
r <- sample(Jester5k, 1000) 
r


###################################################
### code chunk number 14: recommenderlab.Rnw:1187-1190
###################################################
rowCounts(r[1,])
as(r[1,], "list")
rowMeans(r[1,])


###################################################
### code chunk number 15: hist1
###################################################
hist(getRatings(r), breaks=100)


###################################################
### code chunk number 16: hist2
###################################################
hist(getRatings(normalize(r)), breaks=100)


###################################################
### code chunk number 17: hist3
###################################################
hist(getRatings(normalize(r, method="Z-score")), breaks=100)


###################################################
### code chunk number 18: hist4
###################################################
hist(rowCounts(r), breaks=50)


###################################################
### code chunk number 19: hist5
###################################################
hist(colMeans(r), breaks=20)


###################################################
### code chunk number 20: recommenderlab.Rnw:1277-1278
###################################################
recommenderRegistry$get_entries(dataType = "realRatingMatrix")


###################################################
### code chunk number 21: recommenderlab.Rnw:1286-1288
###################################################
r <- Recommender(Jester5k[1:1000], method = "POPULAR")
r


###################################################
### code chunk number 22: recommenderlab.Rnw:1292-1294
###################################################
names(getModel(r))
getModel(r)$topN


###################################################
### code chunk number 23: recommenderlab.Rnw:1309-1311
###################################################
recom <- predict(r, Jester5k[1001:1002], n=5)
recom


###################################################
### code chunk number 24: recommenderlab.Rnw:1316-1317
###################################################
as(recom, "list")


###################################################
### code chunk number 25: recommenderlab.Rnw:1323-1326
###################################################
recom3 <- bestN(recom, n = 3)
recom3
as(recom3, "list")


###################################################
### code chunk number 26: recommenderlab.Rnw:1334-1337
###################################################
recom <- predict(r, Jester5k[1001:1002], type="ratings")
recom
as(recom, "matrix")[,1:10]


###################################################
### code chunk number 27: recommenderlab.Rnw:1357-1360
###################################################
e <- evaluationScheme(Jester5k[1:1000], method="split", train=0.9, 
    given=15, goodRating=5)
e


###################################################
### code chunk number 28: recommenderlab.Rnw:1366-1371
###################################################
r1 <- Recommender(getData(e, "train"), "UBCF")
r1

r2 <- Recommender(getData(e, "train"), "IBCF")
r2


###################################################
### code chunk number 29: recommenderlab.Rnw:1378-1382
###################################################
p1 <- predict(r1, getData(e, "known"), type="ratings")
p1
p2 <- predict(r2, getData(e, "known"), type="ratings")
p2


###################################################
### code chunk number 30: recommenderlab.Rnw:1388-1394
###################################################
error <- rbind(
  calcPredictionAccuracy(p1, getData(e, "unknown")),
  calcPredictionAccuracy(p2, getData(e, "unknown"))
)
rownames(error) <- c("UBCF","IBCF")
error


###################################################
### code chunk number 31: recommenderlab.Rnw:1407-1410
###################################################
scheme <- evaluationScheme(Jester5k[1:1000], method="cross", k=4, given=3,
    goodRating=5)
scheme


###################################################
### code chunk number 32: recommenderlab.Rnw:1417-1420
###################################################
results <- evaluate(scheme, method="POPULAR", type = "topNList", 
  n=c(1,3,5,10,15,20))
results


###################################################
### code chunk number 33: recommenderlab.Rnw:1431-1432
###################################################
getConfusionMatrix(results)[[1]]


###################################################
### code chunk number 34: recommenderlab.Rnw:1444-1445
###################################################
avg(results)


###################################################
### code chunk number 35: roc1
###################################################
plot(results, annotate=TRUE)


###################################################
### code chunk number 36: precrec1
###################################################
plot(results, "prec/rec", annotate=TRUE)


###################################################
### code chunk number 37: recommenderlab.Rnw:1491-1507
###################################################
set.seed(2016)
scheme <- evaluationScheme(Jester5k[1:1000], method="split", train = .9, 
  k=1, given=-5, goodRating=5)
scheme

algorithms <- list(
  "random items" = list(name="RANDOM", param=NULL),
  "popular items" = list(name="POPULAR", param=NULL),
  "user-based CF" = list(name="UBCF", param=list(nn=50)),
  "item-based CF" = list(name="IBCF", param=list(k=50)),
  "SVD approximation" = list(name="SVD", param=list(approxRank = 50))
)

## run algorithms
results <- evaluate(scheme, algorithms, type = "topNList", 
  n=c(1, 3, 5, 10, 15, 20))


###################################################
### code chunk number 38: recommenderlab.Rnw:1512-1513
###################################################
results


###################################################
### code chunk number 39: recommenderlab.Rnw:1519-1521
###################################################
names(results)
results[["user-based CF"]]


###################################################
### code chunk number 40: roc2
###################################################
plot(results, annotate=c(1,3), legend="bottomright")


###################################################
### code chunk number 41: precrec2
###################################################
plot(results, "prec/rec", annotate=3, legend="topleft")


###################################################
### code chunk number 42: recommenderlab.Rnw:1564-1566
###################################################
## run algorithms
results <- evaluate(scheme, algorithms, type = "ratings") 


###################################################
### code chunk number 43: recommenderlab.Rnw:1571-1572
###################################################
results


###################################################
### code chunk number 44: real
###################################################
plot(results, ylim = c(0,100))


###################################################
### code chunk number 45: recommenderlab.Rnw:1595-1601
###################################################
Jester_binary <- binarize(Jester5k, minRating=5)
Jester_binary <- Jester_binary[rowCounts(Jester_binary)>20]
Jester_binary
scheme_binary <- evaluationScheme(Jester_binary[1:1000], 
	method="split", train=.9, k=1, given=3)
scheme_binary


###################################################
### code chunk number 46: recommenderlab.Rnw:1604-1606
###################################################
results_binary <- evaluate(scheme_binary, algorithms, 
  type = "topNList", n=c(1,3,5,10,15,20))


###################################################
### code chunk number 47: roc3
###################################################
plot(results_binary, annotate=c(1,3), legend="bottomright")


