### R code from vignette source 'ctree.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 70, SweaveHooks = list(leftpar = 
    function() par(mai = par("mai") * c(1, 1.1, 1, 1))))
require("partykit")
require("coin")
require("strucchange")
set.seed(290875)


###################################################
### code chunk number 2: party-max
###################################################
ctree_control(teststat = "max")


###################################################
### code chunk number 3: party-max
###################################################
ctree_control(teststat = "quad")


###################################################
### code chunk number 4: party-Bonf
###################################################
ctree_control(testtype = "Bonferroni")


###################################################
### code chunk number 5: party-minsplit
###################################################
ctree_control(minsplit = 20)


###################################################
### code chunk number 6: party-maxsurr
###################################################
ctree_control(maxsurrogate = 3)


###################################################
### code chunk number 7: party-data
###################################################
ls <- data.frame(y = gl(3, 50, labels = c("A", "B", "C")), 
                 x1 = rnorm(150) + rep(c(1, 0, 0), c(50, 50, 50)), 
                 x2 = runif(150))


###################################################
### code chunk number 8: party-formula
###################################################
library("partykit")
ctree(y ~ x1 + x2, data = ls)


###################################################
### code chunk number 9: party-fitted
###################################################
ct <- ctree(y ~ x1 + x2, data = ls)


###################################################
### code chunk number 10: party-print
###################################################
ct


###################################################
### code chunk number 11: party-plot
###################################################
plot(ct)


###################################################
### code chunk number 12: party-nodes
###################################################
ct[1]


###################################################
### code chunk number 13: party-nodelist
###################################################
class(ct[1])


###################################################
### code chunk number 14: party-predict
###################################################
predict(ct, newdata = ls)


###################################################
### code chunk number 15: party-treeresponse
###################################################
predict(ct, newdata = ls[c(1, 51, 101),], type = "prob")


###################################################
### code chunk number 16: party-where
###################################################
predict(ct, newdata = ls[c(1,51,101),], type = "node")


###################################################
### code chunk number 17: party-sctest
###################################################
library("strucchange")
sctest(ct)


###################################################
### code chunk number 18: treepipit-ctree
###################################################
data("treepipit", package = "coin")
tptree <- ctree(counts ~ ., data = treepipit)


###################################################
### code chunk number 19: treepipit-plot
###################################################
plot(tptree, terminal_panel = node_barplot)


###################################################
### code chunk number 20: treepipit-x
###################################################
p <- info_node(node_party(tptree))$p.value
n <- table(predict(tptree, type = "node"))


###################################################
### code chunk number 21: glaucoma-ctree
###################################################
data("GlaucomaM", package = "TH.data")
gtree <- ctree(Class ~ ., data = GlaucomaM)


###################################################
### code chunk number 22: glaucoma-x
###################################################
sp <- split_node(node_party(gtree))$varID


###################################################
### code chunk number 23: glaucoma-plot
###################################################
plot(gtree)


###################################################
### code chunk number 24: glaucoma-plot-inner
###################################################
plot(gtree, inner_panel = node_barplot, 
     edge_panel = function(...) invisible(), tnex = 1)


###################################################
### code chunk number 25: glaucoma-plot2 (eval = FALSE)
###################################################
## plot(gtree)


###################################################
### code chunk number 26: glaucoma-plot-inner (eval = FALSE)
###################################################
## plot(gtree, inner_panel = node_barplot, 
##      edge_panel = function(...) invisible(), tnex = 1)


###################################################
### code chunk number 27: glaucoma-prediction
###################################################
table(predict(gtree), GlaucomaM$Class)


###################################################
### code chunk number 28: glaucoma-classprob
###################################################
prob <- predict(gtree, type = "prob")[,1] + 
                runif(nrow(GlaucomaM), min = -0.01, max = 0.01)
splitvar <- character_split(split_node(node_party(gtree)), 
                            data = data_party(gtree))$name
plot(GlaucomaM[[splitvar]], prob, 
     pch = as.numeric(GlaucomaM$Class), ylab = "Conditional Class Prob.",
     xlab = splitvar)
abline(v = split_node(node_party(gtree))$breaks, lty = 2)
legend(0.15, 0.7, pch = 1:2, legend = levels(GlaucomaM$Class), bty = "n")


###################################################
### code chunk number 29: GBSGS-ctree
###################################################
data("GBSG2", package = "TH.data")  
library("survival")
(stree <- ctree(Surv(time, cens) ~ ., data = GBSG2))


###################################################
### code chunk number 30: GBSG2-plot
###################################################
plot(stree)


###################################################
### code chunk number 31: GBSG2-KM
###################################################
pn <- predict(stree, newdata = GBSG2[1:2,], type = "node")
n <- predict(stree, type = "node")
survfit(Surv(time, cens) ~ 1, data = GBSG2, subset = (n == pn[1]))
survfit(Surv(time, cens) ~ 1, data = GBSG2, subset = (n == pn[2]))


###################################################
### code chunk number 32: mammo-ctree
###################################################
data("mammoexp", package = "TH.data")
mtree <- ctree(ME ~ ., data = mammoexp)


###################################################
### code chunk number 33: mammo-plot
###################################################
plot(mtree)


###################################################
### code chunk number 34: spider-ctree
###################################################
data("HuntingSpiders", package = "partykit")
sptree <- ctree(arct.lute + pard.lugu + zora.spin + pard.nigr +
  pard.pull + aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce +
  alop.fabr + arct.peri ~ herbs + reft + moss + sand + twigs + water,
  data = HuntingSpiders, teststat = "max", minsplit = 5)


###################################################
### code chunk number 35: spider-plot1
###################################################
plot(sptree, terminal_panel = node_barplot)


###################################################
### code chunk number 36: spider-plot2
###################################################
plot(sptree)


