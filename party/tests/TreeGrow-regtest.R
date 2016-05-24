
set.seed(290875)
library("party")
if (!require("TH.data"))
    stop("cannot load package TH.data")
if (!require("coin"))
    stop("cannot load package coin")

### get rid of the NAMESPACE
attach(asNamespace("party"))

gtctrl <- new("GlobalTestControl")
tlev <- levels(gtctrl@testtype)

data(GlaucomaM, package = "TH.data")
gtree <- ctree(Class ~ ., data = GlaucomaM)
tree <- gtree@tree
stopifnot(isequal(tree[[5]][[3]], 0.059))
predict(gtree)

# print(tree)

stump <- ctree(Class ~ ., data = GlaucomaM, 
               control = ctree_control(stump = TRUE))
print(stump)

data(treepipit, package = "coin")

tr <- ctree(counts ~ ., data = treepipit)
tr
plot(tr)


data(GlaucomaM, package = "TH.data")

tr <- ctree(Class ~ ., data = GlaucomaM)
tr
plot(tr)

data(GBSG2, package = "TH.data")  

GBSG2tree <- ctree(Surv(time, cens) ~ ., data = GBSG2)
GBSG2tree
plot(GBSG2tree)
plot(GBSG2tree, terminal_panel = node_surv(GBSG2tree))
survfit(Surv(time, cens) ~ as.factor(GBSG2tree@where), data = GBSG2)
names(GBSG2)

tr <- ctree(Surv(time, cens) ~ ., data = GBSG2, 
            control = ctree_control(teststat = "max", 
                                    testtype = "Univariate"))
tr
plot(tr)

data("mammoexp", package = "TH.data")
attr(mammoexp$ME, "scores") <- 1:3   
attr(mammoexp$SYMPT, "scores") <- 1:4
attr(mammoexp$DECT, "scores") <- 1:3 
names(mammoexp)[names(mammoexp) == "SYMPT"] <- "symptoms"
names(mammoexp)[names(mammoexp) == "PB"] <- "benefit"

names(mammoexp)
tr <- ctree(ME ~ ., data = mammoexp)
tr
plot(tr)

treeresponse(tr, newdata = mammoexp[1:5,])

### check different user interfaces
data("iris")
x <- as.matrix(iris[,colnames(iris) != "Species"])
y <- iris[,"Species"]
newx <- x

ls <- LearningSample(x, y)
p1 <- unlist(treeresponse(ctree(Species ~ ., data = iris), newdata = as.data.frame(newx)))
p2 <- unlist(treeresponse(ctreefit(ls, control = ctree_control()), newdata = as.matrix(newx)))
stopifnot(identical(max(abs(p1 - p2)), 0))

set.seed(29)
p1 <- unlist(treeresponse(cforestfit(ls, control = cforest_unbiased(mtry = 1)), newdata = as.matrix(newx)))
set.seed(29)
p2 <- unlist(treeresponse(cforest(Species ~ ., data = iris, control = cforest_unbiased(mtry = 1)), 
             newdata = as.data.frame(newx)))
stopifnot(identical(max(abs(p1 - p2)), 0))
