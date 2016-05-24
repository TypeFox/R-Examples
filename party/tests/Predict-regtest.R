
set.seed(290875)
library("party")
if (!require("TH.data"))
    stop("cannot load package TH.data")
if (!require("coin"))
    stop("cannot load package coin")

### get rid of the NAMESPACE
attach(asNamespace("party"))

data(treepipit, package = "coin")
ct <- ctree(counts ~ ., data = treepipit)
stopifnot(isequal(predict(ct), predict(ct, newdata = treepipit)))


data(GlaucomaM, package = "TH.data")
ct <- ctree(Class ~ ., data = GlaucomaM)
stopifnot(isequal(predict(ct), predict(ct, newdata = GlaucomaM)))
stopifnot(isequal(predict(ct, type = "prob"), predict(ct, type = "prob", 
                  newdata = GlaucomaM)))
stopifnot(isequal(predict(ct, type = "node"), predict(ct, type = "node", 
                  newdata = GlaucomaM)))
stopifnot(isequal(predict(ct, type = "prob"), treeresponse(ct)))

data("GBSG2", package = "TH.data")  

GBSG2tree <- ctree(Surv(time, cens) ~ ., data = GBSG2)
stopifnot(isequal(GBSG2tree@predict_response(), 
          GBSG2tree@predict_response(newdata = GBSG2)))
stopifnot(isequal(GBSG2tree@cond_distr_response(), 
          GBSG2tree@cond_distr_response(newdata = GBSG2)))

data("mammoexp", package = "TH.data")
attr(mammoexp$ME, "scores") <- 1:3   
attr(mammoexp$SYMPT, "scores") <- 1:4
attr(mammoexp$DECT, "scores") <- 1:3 
names(mammoexp)[names(mammoexp) == "SYMPT"] <- "symptoms"
names(mammoexp)[names(mammoexp) == "PB"] <- "benefit"

names(mammoexp)
mtree <- ctree(ME ~ ., data = mammoexp)
stopifnot(isequal(predict(mtree), predict(mtree, newdata = mammoexp)))
stopifnot(isequal(predict(mtree), predict(mtree, newdata = mammoexp)))
