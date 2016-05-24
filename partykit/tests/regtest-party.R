## load package and fix seed
library("partykit")
set.seed(1)

## rpart: kyphosis data
library("rpart")
data("kyphosis", package = "rpart")
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
pfit <- as.party(fit)
all(predict(pfit, newdata = kyphosis, type = "node") == fit$where)

## J48: iris data
library("RWeka")
data("iris", package = "datasets")
itree <- J48(Species ~ ., data = iris)
pitree <- as.party(itree)
stopifnot(all(predict(pitree) == predict(pitree, newdata = iris[, 3:4])))

all.equal(predict(itree, type = "prob", newdata = iris),
          predict(pitree, type = "prob", newdata = iris))
all.equal(predict(itree,  newdata = iris),
          predict(pitree, newdata = iris))

## rpart/J48: GlaucomaM data
data("GlaucomaM", package = "TH.data")
w <- runif(nrow(GlaucomaM))
fit <- rpart(Class ~ ., data = GlaucomaM, weights = w)
pfit <- as.party(fit)
all(predict(pfit, type = "node") == fit$where)
tmp <- GlaucomaM[sample(1:nrow(GlaucomaM), 100),]
all.equal(predict(fit, type = "prob", newdata = tmp), predict(pfit, type = "prob", newdata = tmp))
all.equal(predict(fit, type = "class", newdata = tmp), predict(pfit, newdata = tmp))
itree <- J48(Class ~ ., data = GlaucomaM)
pitree <- as.party(itree)
all.equal(predict(itree, newdata = tmp, type = "prob"), predict(pitree, newdata = tmp, type = "prob"))

## rpart: airquality data
data("airquality")
aq <- subset(airquality, !is.na(Ozone))
w <- runif(nrow(aq), max = 3)
aqr <- rpart(Ozone ~ ., data = aq, weights = w)
aqp <- as.party(aqr)
tmp <- subset(airquality, is.na(Ozone))
all.equal(predict(aqr, newdata = tmp), predict(aqp, newdata = tmp))

## rpart: GBSG2 data
data("GBSG2", package = "TH.data")
library("survival")
fit <- rpart(Surv(time, cens) ~ ., data = GBSG2)
pfit <- as.party(fit)
pfit$fitted
predict(pfit, newdata = GBSG2[1:100,], type = "prob")
predict(pfit, newdata = GBSG2[1:100,], type = "response")

### multiple responses
f <- fitted(pfit)
f[["(response)"]] <- data.frame(srv = f[["(response)"]], hansi = runif(nrow(f)))
mp <- party(node_party(pfit), fitted = f, data = pfit$data)
class(mp) <- c("constparty", "party")
predict(mp, newdata = GBSG2[1:10,])

### pruning
## create party
data("WeatherPlay", package = "partykit")
py <- party(
  partynode(1L,
    split = partysplit(1L, index = 1:3),
    kids = list(
      partynode(2L,
        split = partysplit(3L, breaks = 75),
        kids = list(
          partynode(3L, info = "yes"),
          partynode(4L, info = "no"))),
      partynode(5L, 
        split = partysplit(3L, breaks = 20),
        kids = list(
          partynode(6L, info = "no"),
          partynode(7L, info = "yes"))),
      partynode(8L,
        split = partysplit(4L, index = 1:2),
        kids = list(
          partynode(9L, info = "yes"),
          partynode(10L, info = "no"))))),
  WeatherPlay)
names(py) <- LETTERS[nodeids(py)]

## print
print(py)
(py5 <- nodeprune(py, 5))
nodeids(py5)
(pyH <- nodeprune(py5, "H"))
nodeids(pyH)
	
ct <- ctree(Species ~ ., data = iris)
nt <- node_party(ctree(Species ~ ., data = iris))
(ctp <- nodeprune(ct, 4)) ### party method
(ntp <- nodeprune(nt, 4)) ### partynode method

### check if both methods do the same
p1 <- predict(party(ntp, 
    data = model.frame(terms(ct), data = iris)), type = "node")
p2 <- predict(ctp, type = "node")
stopifnot(max(abs(p1 - p2)) == 0)

names(ct) <- LETTERS[nodeids(ct)]
(ctp <- nodeprune(ct, "D"))

table(predict(ct, type = "node"), 
      predict(ctp, type = "node"))

(ct <- nodeprune(ct, names(ct)[names(ct) != "A"]))
table(predict(ct, type = "node"))

nodeprune(ct, "B")
nodeprune(ct, "C")
nodeprune(ct, "A")

### check different predict flavours for numeric responses
x <- runif(100)
dd <- data.frame(y = rnorm(length(x), mean = 2 * (x < .5)), x = x)
ct <- ctree(y ~ x, data = dd)
nd <- data.frame(x = (1:9) / 10)
predict(ct, newdata = nd, type = "node")
predict(ct, newdata = nd, type = "response")
predict(ct, newdata = nd, type = "prob")
predict(ct, newdata = nd, type = "quantile")
predict(ct, newdata = nd, type = "quantile", at = NULL)
predict(ct, newdata = nd, type = "density")
predict(ct, newdata = nd, type = "density", at = (1:9) / 10)
