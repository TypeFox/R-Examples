### R code from vignette source 'constparty.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 70)
library("partykit")
set.seed(290875)


###################################################
### code chunk number 2: Titanic
###################################################
data("Titanic", package = "datasets")
ttnc <- as.data.frame(Titanic)
ttnc <- ttnc[rep(1:nrow(ttnc), ttnc$Freq), 1:4]
names(ttnc)[2] <- "Gender"


###################################################
### code chunk number 3: rpart
###################################################
library("rpart")
(rp <- rpart(Survived ~ ., data = ttnc))


###################################################
### code chunk number 4: rpart-party
###################################################
(party_rp <- as.party(rp))


###################################################
### code chunk number 5: rpart-plot-orig
###################################################
plot(rp)
text(rp)


###################################################
### code chunk number 6: rpart-plot
###################################################
plot(party_rp)


###################################################
### code chunk number 7: rpart-pred
###################################################
all.equal(predict(rp), predict(party_rp, type = "prob"), 
  check.attributes = FALSE)


###################################################
### code chunk number 8: rpart-fitted
###################################################
str(fitted(party_rp))


###################################################
### code chunk number 9: rpart-prob
###################################################
prop.table(do.call("table", fitted(party_rp)), 1)


###################################################
### code chunk number 10: J48
###################################################
library("RWeka")
(j48 <- J48(Survived ~ ., data = ttnc))


###################################################
### code chunk number 11: J48-party
###################################################
(party_j48 <- as.party(j48))


###################################################
### code chunk number 12: J48-plot
###################################################
plot(party_j48)


###################################################
### code chunk number 13: J48-pred
###################################################
all.equal(predict(j48, type = "prob"), predict(party_j48, type = "prob"),
  check.attributes = FALSE)


###################################################
### code chunk number 14: PMML-Titantic
###################################################
ttnc_pmml <- file.path(system.file("pmml", package = "partykit"),
  "ttnc.pmml")
(ttnc_quest <- pmmlTreeModel(ttnc_pmml))


###################################################
### code chunk number 15: PMML-Titanic-plot1
###################################################
plot(ttnc_quest)


###################################################
### code chunk number 16: ttnc2-reorder
###################################################
ttnc2 <- ttnc[, names(ttnc_quest$data)]
for(n in names(ttnc2)) {
  if(is.factor(ttnc2[[n]])) ttnc2[[n]] <- factor(
    ttnc2[[n]], levels = levels(ttnc_quest$data[[n]]))
}


###################################################
### code chunk number 17: PMML-Titanic-augmentation
###################################################
ttnc_quest2 <- party(ttnc_quest$node,
  data = ttnc2,
  fitted = data.frame(
    "(fitted)" = predict(ttnc_quest, ttnc2, type = "node"),
    "(response)" = ttnc2$Survived,
    check.names = FALSE),
  terms = terms(Survived ~ ., data = ttnc2)
)
ttnc_quest2 <- as.constparty(ttnc_quest2)


###################################################
### code chunk number 18: PMML-Titanic-plot2
###################################################
plot(ttnc_quest2)


###################################################
### code chunk number 19: PMML-write
###################################################
library("pmml")
tfile <- tempfile()
write(toString(pmml(rp)), file = tfile)


###################################################
### code chunk number 20: PMML-read
###################################################
(party_pmml <- pmmlTreeModel(tfile))
all.equal(predict(party_rp, newdata = ttnc, type = "prob"), 
  predict(party_pmml, newdata = ttnc, type = "prob"),
  check.attributes = FALSE)


###################################################
### code chunk number 21: mytree-1
###################################################
findsplit <- function(response, data, weights, alpha = 0.01) {

  ## extract response values from data
  y <- factor(rep(data[[response]], weights))

  ## perform chi-squared test of y vs. x
  mychisqtest <- function(x) {
    x <- factor(x)
    if(length(levels(x)) < 2) return(NA)
    ct <- suppressWarnings(chisq.test(table(y, x), correct = FALSE))
    pchisq(ct$statistic, ct$parameter, log = TRUE, lower.tail = FALSE)
  }
  xselect <- which(names(data) != response)
  logp <- sapply(xselect, function(i) mychisqtest(rep(data[[i]], weights)))
  names(logp) <- names(data)[xselect]

  ## Bonferroni-adjusted p-value small enough?
  if(all(is.na(logp))) return(NULL)
  minp <- exp(min(logp, na.rm = TRUE))
  minp <- 1 - (1 - minp)^sum(!is.na(logp))
  if(minp > alpha) return(NULL)

  ## for selected variable, search for split minimizing p-value  
  xselect <- xselect[which.min(logp)]
  x <- rep(data[[xselect]], weights)

  ## set up all possible splits in two kid nodes
  lev <- levels(x[drop = TRUE])
  if(length(lev) == 2) {
    splitpoint <- lev[1]
  } else {
    comb <- do.call("c", lapply(1:(length(lev) - 2),
      function(x) combn(lev, x, simplify = FALSE)))
    xlogp <- sapply(comb, function(q) mychisqtest(x %in% q))
    splitpoint <- comb[[which.min(xlogp)]]
  }

  ## split into two groups (setting groups that do not occur to NA)
  splitindex <- !(levels(data[[xselect]]) %in% splitpoint)
  splitindex[!(levels(data[[xselect]]) %in% lev)] <- NA_integer_
  splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L

  ## return split as partysplit object
  return(partysplit(varid = as.integer(xselect),
    index = splitindex,
    info = list(p.value = 1 - (1 - exp(logp))^sum(!is.na(logp)))))
}


###################################################
### code chunk number 22: mytree-2
###################################################
growtree <- function(id = 1L, response, data, weights, minbucket = 30) {

  ## for less than 30 observations stop here
  if (sum(weights) < minbucket) return(partynode(id = id))

  ## find best split
  sp <- findsplit(response, data, weights)
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id))

  ## actually split the data
  kidids <- kidids_split(sp, data = data)

  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1:length(kids)) {
  ## select observations for current node
  w <- weights
  w[kidids != kidid] <- 0
  ## get next node id
  if (kidid > 1) {
    myid <- max(nodeids(kids[[kidid - 1]]))
  } else {
    myid <- id
  }
  ## start recursion on this daugther node
  kids[[kidid]] <- growtree(id = as.integer(myid + 1), response, data, w)
  }

  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
    info = list(p.value = min(info_split(sp)$p.value, na.rm = TRUE))))
}


###################################################
### code chunk number 23: mytree-3
###################################################
mytree <- function(formula, data, weights = NULL) {

  ## name of the response variable
  response <- all.vars(formula)[1]
  ## data without missing values, response comes last
  data <- data[complete.cases(data), c(all.vars(formula)[-1], response)]
  ## data is factors only
  stopifnot(all(sapply(data, is.factor)))

  if (is.null(weights)) weights <- rep(1L, nrow(data))
  ## weights are case weights, i.e., integers
  stopifnot(length(weights) == nrow(data) &
    max(abs(weights - floor(weights))) < .Machine$double.eps)

  ## grow tree
  nodes <- growtree(id = 1L, response, data, weights)

  ## compute terminal node number for each observation
  fitted <- fitted_node(nodes, data = data)
  ## return rich constparty object
  ret <- party(nodes, data = data,
    fitted = data.frame("(fitted)" = fitted,
                        "(response)" = data[[response]],
                        "(weights)" = weights,
                        check.names = FALSE),
    terms = terms(formula))
  as.constparty(ret)
}


###################################################
### code chunk number 24: mytree-4
###################################################
(myttnc <- mytree(Survived ~ Class + Age + Gender, data = ttnc))


###################################################
### code chunk number 25: mytree-5
###################################################
plot(myttnc)


###################################################
### code chunk number 26: mytree-pval
###################################################
nid <- nodeids(myttnc)
iid <- nid[!(nid %in% nodeids(myttnc, terminal = TRUE))]
(pval <- unlist(nodeapply(myttnc, ids = iid,
  FUN = function(n) info_node(n)$p.value)))


###################################################
### code chunk number 27: mytree-nodeprune
###################################################
myttnc2 <- nodeprune(myttnc, ids = iid[pval > 1e-5])


###################################################
### code chunk number 28: mytree-nodeprune-plot
###################################################
plot(myttnc2)


###################################################
### code chunk number 29: mytree-glm
###################################################
logLik(glm(Survived ~ Class + Age + Gender, data = ttnc, 
           family = binomial()))


###################################################
### code chunk number 30: mytree-bs
###################################################
bs <- rmultinom(25, nrow(ttnc), rep(1, nrow(ttnc)) / nrow(ttnc))


###################################################
### code chunk number 31: mytree-ll
###################################################
bloglik <- function(prob, weights)
    sum(weights * dbinom(ttnc$Survived == "Yes", size = 1, 
                         prob[,"Yes"], log = TRUE))


###################################################
### code chunk number 32: mytree-bsll
###################################################
f <- function(w) {
    tr <- mytree(Survived ~ Class + Age + Gender, data = ttnc, weights = w)
    bloglik(predict(tr, newdata = ttnc, type = "prob"), as.numeric(w == 0))
}
apply(bs, 2, f)


###################################################
### code chunk number 33: mytree-node
###################################################
nttnc <- expand.grid(Class = levels(ttnc$Class),
  Gender = levels(ttnc$Gender), Age = levels(ttnc$Age))
nttnc


###################################################
### code chunk number 34: mytree-prob
###################################################
predict(myttnc, newdata = nttnc, type = "node")
predict(myttnc, newdata = nttnc, type = "response")
predict(myttnc, newdata = nttnc, type = "prob")


###################################################
### code chunk number 35: mytree-FUN
###################################################
predict(myttnc, newdata = nttnc, FUN = function(y, w)
  rank(table(rep(y, w))))


