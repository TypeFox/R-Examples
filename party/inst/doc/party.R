### R code from vignette source 'party.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 70, SweaveHooks = list(leftpar = 
    function() par(mai = par("mai") * c(1, 1.1, 1, 1))))
require("party")
require("coin")
set.seed(290875)


###################################################
### code chunk number 2: party-data
###################################################
ls <- data.frame(y = gl(3, 50, labels = c("A", "B", "C")), x1 = rnorm(150) + rep(c(1, 0, 0), c(50, 50,
50)), x2 = runif(150))


###################################################
### code chunk number 3: party-formula
###################################################
library("party")
ctree(y ~ x1 + x2, data = ls)


###################################################
### code chunk number 4: party-xtrafo
###################################################
ctree(y ~ x1 + x2, data = ls, xtrafo = function(data) trafo(data,
numeric_trafo = rank))


###################################################
### code chunk number 5: party-Bonf
###################################################
ctree_control(testtype = "Bonferroni")


###################################################
### code chunk number 6: party-MC
###################################################
ctree_control(testtype = "MonteCarlo")           


###################################################
### code chunk number 7: party-ss
###################################################
ctree_control(savesplitstats = TRUE)


###################################################
### code chunk number 8: party-minsplit
###################################################
ctree_control(minsplit = 20)


###################################################
### code chunk number 9: party-maxsurr
###################################################
ctree_control(maxsurrogate = 3)


###################################################
### code chunk number 10: party-fitted
###################################################
ct <- ctree(y ~ x1 + x2, data = ls)


###################################################
### code chunk number 11: party-print
###################################################
ct


###################################################
### code chunk number 12: party-plot
###################################################
plot(ct)


###################################################
### code chunk number 13: party-nodes
###################################################
nodes(ct, 1)


###################################################
### code chunk number 14: party-nodelist
###################################################
names(nodes(ct, 1)[[1]])


###################################################
### code chunk number 15: party-Predict
###################################################
Predict(ct, newdata = ls)


###################################################
### code chunk number 16: party-treeresponse
###################################################
treeresponse(ct, newdata = ls[c(1, 51, 101),])


###################################################
### code chunk number 17: party-where
###################################################
where(ct, newdata = ls[c(1,51,101),])


###################################################
### code chunk number 18: treepipit-ctree
###################################################
data("treepipit", package = "coin")
tptree <- ctree(counts ~ ., data = treepipit)


###################################################
### code chunk number 19: treepipit-plot
###################################################
plot(tptree, terminal_panel = node_hist(tptree, breaks = 0:6-0.5, ymax = 65, 
horizontal = FALSE, freq = TRUE))


###################################################
### code chunk number 20: treepipit-x
###################################################
x <- tptree@tree


###################################################
### code chunk number 21: glaucoma-ctree
###################################################
data("GlaucomaM", package = "TH.data")
gtree <- ctree(Class ~ ., data = GlaucomaM)


###################################################
### code chunk number 22: glaucoma-x
###################################################
x <- gtree@tree


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
### code chunk number 27: glaucoma-split-plot
###################################################
getOption("SweaveHooks")[["leftpar"]]()
cex <- 1.6
inner <- nodes(gtree, c(1, 2, 5))
layout(matrix(1:length(inner), ncol = length(inner)))
out <- sapply(inner, function(i) {
    splitstat <- i$psplit$splitstatistic
    x <- GlaucomaM[[i$psplit$variableName]][splitstat > 0]
    plot(x, splitstat[splitstat > 0], main = paste("Node", i$nodeID),
         xlab = i$psplit$variableName, ylab = "Statistic", ylim = c(0, 10),
         cex.axis = cex, cex.lab = cex, cex.main = cex)
    abline(v = i$psplit$splitpoint, lty = 3)
})


###################################################
### code chunk number 28: glaucoma-prediction
###################################################
table(Predict(gtree), GlaucomaM$Class)


###################################################
### code chunk number 29: glaucoma-classprob
###################################################
prob <- sapply(treeresponse(gtree), function(x) x[1]) +
               runif(nrow(GlaucomaM), min = -0.01, max = 0.01)
splitvar <- nodes(gtree, 1)[[1]]$psplit$variableName
plot(GlaucomaM[[splitvar]], prob, 
     pch = as.numeric(GlaucomaM$Class), ylab = "Conditional Class Prob.",
     xlab = splitvar)
abline(v = nodes(gtree, 1)[[1]]$psplit$splitpoint, lty = 2)
legend(0.15, 0.7, pch = 1:2, legend = levels(GlaucomaM$Class), bty = "n")


###################################################
### code chunk number 30: GBSGS-ctree
###################################################
data("GBSG2", package = "TH.data")  
stree <- ctree(Surv(time, cens) ~ ., data = GBSG2)


###################################################
### code chunk number 31: GBSG2-plot
###################################################
plot(stree)


###################################################
### code chunk number 32: GBSG2-KM
###################################################
treeresponse(stree, newdata = GBSG2[1:2,])


###################################################
### code chunk number 33: mammo-ctree
###################################################
data("mammoexp", package = "TH.data")
mtree <- ctree(ME ~ ., data = mammoexp)


###################################################
### code chunk number 34: mammo-plot
###################################################
plot(mtree)


###################################################
### code chunk number 35: spider-ctree (eval = FALSE)
###################################################
## data("spider", package = "mvpart")   
## sptree <- ctree(arct.lute + pard.lugu + zora.spin + pard.nigr + pard.pull +
##            aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce + 
##            alop.fabr + arct.peri ~ herbs+reft+moss+sand+twigs+water, 
##            controls = ctree_control(teststat = "max", testtype = "MonteCarlo", nresample = 9999, 
##            minsplit = 5, mincriterion = 0.9), data = spider)
## sptree@tree$criterion
## 
## library("coin")
## data("spider", package = "mvpart")
## it <- independence_test(arct.lute + pard.lugu + zora.spin + pard.nigr + pard.pull +
##                        aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce +   
##                        alop.fabr + arct.peri ~ herbs+reft+moss+sand+twigs+water, 
##                        data = spider, distribution = approximate(B = 19999))
## statistic(it, "standardized")
## pvalue(it)


###################################################
### code chunk number 36: spider-plot (eval = FALSE)
###################################################
## plot(sptree, terminal_panel = node_terminal)


