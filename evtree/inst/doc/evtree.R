### R code from vignette source 'evtree.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("rpart")
library("evtree")
library("lattice")
data("BBBClub", package = "evtree")
cache <- FALSE


###################################################
### code chunk number 2: chess22
###################################################
X1 <- rep(seq(0.25, 1.75, 0.5), each = 4)
X2 <- rep(seq(0.25, 1.75, 0.5), 4)
Y <- rep(1, 16)
Y[(X1 < 1 & X2 < 1) | (X1 > 1 & X2 > 1)] <- 2
Y <- factor(Y, labels = c("O", "X"))
chess22 <- data.frame(Y, X1, X2)
set.seed(1090)
print(evtree(Y ~ ., data = chess22, minbucket = 1, minsplit = 2))


###################################################
### code chunk number 3: chess22-plot
###################################################
par(mar = c(4, 4, 1, 1))
plot(X2 ~ X1, data = chess22, xlim = c(0, 2), ylim = c(0, 2), pch = c(1, 4)[Y], col = c("black", "slategray")[Y])


###################################################
### code chunk number 4: BBBClub-rpart-ctree (eval = FALSE)
###################################################
## data("BBBClub", package = "evtree")
## library("rpart")
## rp  <- as.party(rpart(choice ~ ., data = BBBClub, minbucket = 10))
## rp2 <- as.party(rpart(choice ~ ., data = BBBClub, minbucket = 10,
##   maxdepth = 2))
## ct  <- ctree(choice ~ ., data = BBBClub, minbucket = 10, mincrit = 0.99)
## ct2 <- ctree(choice ~ ., data = BBBClub, minbucket = 10, mincrit = 0.99,
##   maxdepth = 2)
## plot(rp)
## plot(ct)


###################################################
### code chunk number 5: BBBClub-evtree (eval = FALSE)
###################################################
## set.seed(1090)
## ev <- evtree(choice ~ ., data = BBBClub, minbucket = 10, maxdepth = 2)


###################################################
### code chunk number 6: BBBClub-cache
###################################################
if(cache & file.exists("BBBClub-trees.rda")) {
load("BBBClub-trees.rda")
} else {
data("BBBClub", package = "evtree")
library("rpart")
rp  <- as.party(rpart(choice ~ ., data = BBBClub, minbucket = 10))
rp2 <- as.party(rpart(choice ~ ., data = BBBClub, minbucket = 10,
  maxdepth = 2))
ct  <- ctree(choice ~ ., data = BBBClub, minbucket = 10, mincrit = 0.99)
ct2 <- ctree(choice ~ ., data = BBBClub, minbucket = 10, mincrit = 0.99,
  maxdepth = 2)
plot(rp)
plot(ct)
set.seed(1090)
ev <- evtree(choice ~ ., data = BBBClub, minbucket = 10, maxdepth = 2)
if(cache) {
  save(rp, rp2, ct, ct2, ev, file = "BBBClub-trees.rda")
} else {
  if(file.exists("BBBClub-trees.rda")) file.remove("BBBClub-trees.rda")
}
}


###################################################
### code chunk number 7: BBBClub-rpart-plot
###################################################
plot(rp)


###################################################
### code chunk number 8: BBBClub-ctree-plot
###################################################
plot(ct)


###################################################
### code chunk number 9: BBBClub-evtree-display
###################################################
plot(ev)
ev


###################################################
### code chunk number 10: BBBClub-evtree-plot
###################################################
plot(ev)


###################################################
### code chunk number 11: evtree-performance
###################################################
mc <- function(obj) 1 - mean(predict(obj) == BBBClub$choice)
evalfun <- function(obj) 2 * nrow(BBBClub) * mc(obj) +
  width(obj) * log(nrow(BBBClub))
trees <- list("evtree" = ev, "rpart" = rp, "ctree" = ct, "rpart2" = rp2,
  "ctree2" = ct2)
round(sapply(trees, function(obj) c("misclassification" = mc(obj),
  "evaluation function" = evalfun(obj))), digits = 3)


###################################################
### code chunk number 12: evtree-structure
###################################################
ftable(tab <- table(evtree = predict(ev), rpart  = predict(rp),
  ctree  = predict(ct), observed = BBBClub$choice))
sapply(c("evtree", "rpart", "ctree"), function(nam) {
  mt <- margin.table(tab, c(match(nam, names(dimnames(tab))), 4))
  c(abs = as.vector(rowSums(mt))[2],
    rel = round(100 * prop.table(mt, 1)[2, 2], digits = 3))
})


###################################################
### code chunk number 13: benchmark-results
###################################################
## load results
rm(list = ls())
for(i in Sys.glob("results/*.RData")) load(i)
for(i in Sys.glob("results_j48/*.RData")) load(i)

## preprocess for reference evtree
preprocess <- function(d, dname = "datasetname", isclassification = TRUE){
    if(isclassification){
        colnames(d) <- c("evtree", "rpart", "ctree", "J48","evtree", "rpart", "ctree", "J48")
        d[, 1:4] <- 1 - d[ ,1:4]
    }else{
    	colnames(d) <- c("evtree", "rpart", "ctree","evtree", "rpart", "ctree")    	
    }
    d <- as.data.frame(d)
	nAlgorithms = dim(d)[2]/2
    for(i in nAlgorithms:1) d[, i] <- d[, i] / d[, 1] * 100
    if(isclassification)  ## for J48 the total number of nodes is used
         d[, nAlgorithms*2] <- d[, nAlgorithms*2] / (d[, nAlgorithms+1]*2+1) * 100
    else
	d[, nAlgorithms*2] <- d[, nAlgorithms*2] / d[, nAlgorithms+1] * 100
    
    for(i in (nAlgorithms*2-1):(nAlgorithms+1)) d[, i] <- d[, i] / d[, nAlgorithms+1] * 100
    x <- d[, 1:nAlgorithms]
    y <- d[, (nAlgorithms+1):(nAlgorithms*2)]
    rval <- reshape(x, idvar="samp", times=names(x), timevar = "alg",varying= list(names(x)), direction="long")
    names(rval)[2] <- "accuracy"
    rval$complexity <- reshape(y, idvar="samp", times=names(y), timevar = "alg",varying= list(names(y)), direction="long")[,2]
    if(isclassification) 
    	rval$alg <- factor(rval$alg, levels = c("evtree", "ctree", "rpart", "J48"))
    else 
    	rval$alg <- factor(rval$alg, levels = c("evtree", "ctree", "rpart"))
    rval$ds <- dname
    rval
}

## collect results for all datasets
r <- rbind(
preprocess(d = cbind(rglass[,1:3], rglass2[,3], rglass[,4:6], rglass2[,4]), dname = "Glass identification", isclassification = TRUE),
preprocess(d = cbind(rheart[,1:3], rheart2[,3], rheart[,4:6], rheart2[,4]), dname = "Statlog heart", isclassification = TRUE),
preprocess(d = cbind(rionosphere[,1:3], rionosphere2[,3], rionosphere[,4:6], rionosphere2[,4]), dname = "Ionosphere", isclassification = TRUE),
preprocess(d = cbind(rmusk[,1:3], rmusk2[,3], rmusk[,4:6], rmusk2[,4]), dname = "Musk", isclassification = TRUE),
preprocess(d = cbind(rbreastcancer[,1:3], rbreastcancer2[,3], rbreastcancer[,4:6], rbreastcancer2[,4]), dname = "Breast cancer database", isclassification = TRUE),
preprocess(d = cbind(rpima[,1:3], rpima2[,3], rpima[,4:6], rpima2[,4]), dname = "Pima Indians diabetes", isclassification = TRUE),
preprocess(d = cbind(rvowel[,1:3], rvowel2[,3], rvowel[,4:6], rvowel2[,4]), dname = "Vowel", isclassification = TRUE),
preprocess(d = cbind(rcredit[,1:3], rcredit2[,3], rcredit[,4:6], rcredit2[,4]), dname = "Statlog German credit", isclassification = TRUE),
preprocess(d = cbind(rcontraceptive[,1:3], rcontraceptive2[,3], rcontraceptive[,4:6], rcontraceptive2[,4]), dname = "Contraceptive method", isclassification = TRUE),
preprocess(d = cbind(rdna[,1:3], rdna2[,3], rdna[,4:6], rdna2[,4]), dname = "DNA", isclassification = TRUE),
preprocess(d = cbind(rspam[,1:3], rspam2[,3], rspam[,4:6], rspam2[,4]), dname = "Spam", isclassification = TRUE),
preprocess(d = cbind(rmagicgamma[,1:3], rmagicgamma2[,3], rmagicgamma[,4:6], rmagicgamma2[,4]), dname = "Magic gamma telescope", isclassification = TRUE),
preprocess(d = rservo, dname = "Servo", isclassification = FALSE),
preprocess(d = rbostonhousing, dname = "Boston housing", isclassification = FALSE),
preprocess(d = rmel0101, dname = "MEL0101", isclassification = FALSE),
preprocess(d = rhdg0202, dname = "HDG0202", isclassification = FALSE),
preprocess(d = rhdg0502, dname = "HDG0502", isclassification = FALSE)
)

r$ds <- factor(r$ds)
r$samp <- factor(r$samp)
r$dssamp <- r$ds:r$samp

## compute multiple comparisons
library("multcomp")
cstats <- function(alg = "evtree", value = "accuracy", data = r) {
  dlab <- rev(unique(data$ds))
  if(alg == "J48"){ 
  	dlab <- dlab[-c(1:5)] ## J48: skip regression datasets
  }
  k <- length(dlab)  
  mean  <- numeric(k)
  lower <- numeric(k)
  upper <- numeric(k)
  names(data)[names(data) == value] <- "value"
  firstDS <- 1
  for(i in 1:k) {
  	dsub <- subset(data, ds == dlab[i])
	dsub$alg <- factor(dsub$alg)
    mod1 <- lm(value ~ alg, data = dsub)
    pt <- glht(mod1, linfct = mcp(alg = "Dunnett"))
    w <- confint(pt)$confint
    d <- which(levels(dsub$alg) == alg) - 1
    mean[i]  <-  w[d]
    lower[i] <-  w[d + length(levels(dsub$alg))-1]
    upper[i] <-  w[d + (length(levels(dsub$alg))-1)*2]
  }
  rval <- data.frame(mean, lower, upper)
  rownames(rval) <- dlab
  return(rval)
}

acc_rpart <- cstats("rpart", "accuracy")
com_rpart <- cstats("rpart", "complexity")
acc_ctree <- cstats("ctree", "accuracy")
com_ctree <- cstats("ctree", "complexity")
acc_J48 <- cstats("J48", "accuracy") 
com_J48 <- cstats("J48", "complexity") 

## function for visualization
ciplot <- function(x, xlim = NULL, main = "", xlab = "", ylab = TRUE) {
  nam <- rownames(x)
  k <- length(nam)
  plot(x$mean, 1:k, xlim = xlim, axes = FALSE, xlab = "", ylab = "", pch = 19)
  arrows(x$lower, 1:k, x$upper, 1:k, angle = 90, length = 0.05, code = 3)
  if(xlab == "") axis(1, labels = FALSE) else axis(1)
  if(ylab) ylab <- nam
  axis(2, at = 1:k, labels = ylab, las = 1, cex = 0.8)  
  axis(2, at = k + 1.5, labels = main, tick = FALSE, las = 1, outer = TRUE, cex.axis = 1.5, xpd = TRUE)
  mtext(xlab, side = 1, line = 3, xpd = TRUE)
  if (NROW(x) >= 17) abline(h = 5.5)
  abline(v = 0, lty = 2)  
  box()
}

## plot the results if evtree vs. rpart and evtree vs. ctree
par(mfrow = c(2, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(cbind(acc_rpart, acc_ctree))
xlim2 <- range(cbind(com_rpart, com_ctree))

ciplot(acc_rpart, xlim = xlim1, main = "rpart", ylab = TRUE, xlab = "")
ciplot(com_rpart, xlim = xlim2, main = "",      ylab = FALSE, xlab = "")
ciplot(acc_ctree, xlim = xlim1, main = "ctree", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_ctree, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")


## plot the results of evtree vs. J48
par(mfrow = c(1, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(acc_J48)
xlim2 <- range(com_J48)

ciplot(acc_J48, xlim = xlim1, main = "J48", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_J48, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")


###################################################
### code chunk number 14: benchmark-plot
###################################################
par(mfrow = c(2, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(cbind(acc_rpart, acc_ctree))
xlim2 <- range(cbind(com_rpart, com_ctree))

ciplot(acc_rpart, xlim = xlim1, main = "rpart", ylab = TRUE, xlab = "")
ciplot(com_rpart, xlim = xlim2, main = "",      ylab = FALSE, xlab = "")
ciplot(acc_ctree, xlim = xlim1, main = "ctree", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_ctree, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")


###################################################
### code chunk number 15: benchmark-plot2
###################################################
par(mfrow = c(1, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(acc_J48)
xlim2 <- range(com_J48)

ciplot(acc_J48, xlim = xlim1, main = "J48", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_J48, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")


###################################################
### code chunk number 16: chessboard
###################################################
chessboard44 <- function(n = 4000, noisevariables = 6, noise = 0) {
  chess44 <- array(0,c(n,noisevariables+3))
  for(i in 1:(noisevariables+2))
      chess44[,i] <- as.numeric(runif(dim(chess44)[1]))*4

   x <- chess44[,1]
   y <- chess44[,2]
   chess44[, ncol(chess44)] <- 0
   for(k in 1:4)  
      chess44[(x <= k & x > k-1 & y <= k & y > k-1), ncol(chess44)] <- 1
   for(k in 1:2)  
      chess44[(x <= k & x > k-1 & y <= k+2 & y > k+1), ncol(chess44)] <- 1
   for(k in 1:2)  
      chess44[(y <= k & y > k-1 & x <= k+2 & x > k+1), ncol(chess44)] <- 1

   if(noise > 0) {
      flipclasslist <- sample(n, n * (noise / 100), replace = FALSE)

      for(i in 1:length(flipclasslist)){
	  if(chess44[flipclasslist[i], ncol(chess44)] == 1)
	      chess44[flipclasslist[i], ncol(chess44)] = 0
	  else if(chess44[flipclasslist[i], ncol(chess44)] == 0)
	      chess44[flipclasslist[i], ncol(chess44)] = 1
      }
  }

  chess44 <- as.data.frame(chess44)
  chess44[,ncol(chess44)] <- as.factor(chess44[,ncol(chess44)])
  names(chess44) <- c(paste("X", 1:8, sep = ""), "Y")
  chess44
}


###################################################
### code chunk number 17: chessboard-plot
###################################################
chess44 <- chessboard44(2000)
plot(X2 ~ X1, data = chess44, xlim = c(0, 4), ylim = c(0, 4), pch = c(1, 4)[Y], col = c("black", "slategray")[Y])


###################################################
### code chunk number 18: chessboard-table
###################################################
library("xtable")
load("./results/chessboard44_0.RData")
load("./results/chessboard44_5.RData")
load("./results/chessboard44_10.RData")
load("./results_j48/chessboard44_0_j48.RData")
load("./results_j48/chessboard44_5_j48.RData")
load("./results_j48/chessboard44_10_j48.RData")

chesstable_means  <- as.data.frame( rbind(apply(rchessboard44_0,2,mean), apply(rchessboard44_5,2,mean) , apply(rchessboard44_10,2,mean) ) ) 
names(chesstable_means) <-  c("\\code{evtree}", "\\code{rpart}", "\\code{ctree}", "\\code{evtree}", "\\code{rpart}", "\\code{ctree}")
chesstable_means[,1:3] <-  format(chesstable_means[,1:3]*100, digits=1, nsmall=1)
chesstable_means[,4:6] <-  format(chesstable_means[,4:6], digits=1, nsmall=1)

chesstable_sd  <- as.data.frame( rbind(apply(rchessboard44_0,2,sd), apply(rchessboard44_5,2,sd) , apply(rchessboard44_10,2,sd) )) 
names(chesstable_sd) <-  c("\\code{evtree}", "\\code{rpart}", "\\code{ctree}", "\\code{evtree}", "\\code{rpart}", "\\code{ctree}")
chesstable_sd[,1:3] <-  format(chesstable_sd[,1:3]*100, digits=1, nsmall=1)
chesstable_sd[,4:6] <-  format(chesstable_sd[,4:6], digits=1, nsmall=1)

chesstable_means2  <- as.data.frame( cbind(rbind(mean(rchessboard44_02[,3]), mean(rchessboard44_52[,3]) , mean(rchessboard44_102[,3])), rbind(
mean(rchessboard44_02[,4]), mean(rchessboard44_52[,4]) , mean(rchessboard44_102[,4])  )) )
chesstable_means2[,1] <-  format(chesstable_means2[,1]*100, digits=1, nsmall=1)
chesstable_means2[,2] <-  format(chesstable_means2[,2], digits=1, nsmall=1)
names(chesstable_means2) <- c("\\code{J48}", "\\code{J48}") 

chesstable_sd2  <- as.data.frame( cbind(rbind(sd(rchessboard44_02[,3]), sd(rchessboard44_52[,3]) , sd(rchessboard44_102[,3])),rbind(
sd(rchessboard44_02[,4]), sd(rchessboard44_52[,4]) , sd(rchessboard44_102[,4])  ))) 
chesstable_sd2[,1] <-  format(chesstable_sd2[,1]*100, digits=1, nsmall=1)
chesstable_sd2[,2] <-  format(chesstable_sd2[,2], digits=1, nsmall=1)

chesstable_means <- cbind(chesstable_means, chesstable_means2)
chesstable_sd <- cbind(chesstable_sd, chesstable_sd2)
chesstable_means <- chesstable_means[, c(1:3, 7, 4:6, 8)]
chesstable_sd <- chesstable_sd[, c(1:3, 7, 4:6, 8)]

chesstable <- chesstable_means
for(j in 1:ncol(chesstable_means)){
	for(i in 1:nrow(chesstable_means)){
		chesstable[i,j] <- paste(chesstable_means[i,j] ,  "(", chesstable_sd[i,j], ")",  sep="")	
	}
}

chesstable <- cbind(rbind("0\\%","5\\%","10\\%"), chesstable)
colnames(chesstable)[1] <- ""
colnames(chesstable)[6:9] <- colnames(chesstable)[2:5]

print(xtable(chesstable,
caption = "Mean (and standard deviation) of accuracy and complexity for simulated $4 \\times 4$ chessboard examples.",
caption.placement= "bottom",
table.placement="b!",
label= "tab:resultsChessboard"), 
comment = FALSE,
include.rownames = FALSE, allign= "rllllll", hline.after=NULL,
sanitize.text.function = identity,
add.to.row=list(pos=list(-1,-1, 0, 3), command=c(
"\\toprule", 
c("\\multicolumn{1}{l}{Noise} & \\multicolumn{4}{l}{Accuracy}  & \\multicolumn{4}{l}{Complexity}\\\\",
"\\midrule",
"\\bottomrule"
)))
)


###################################################
### code chunk number 19: parameters
###################################################
## load results
for(i in Sys.glob("results_parameter/*.RData")) load(i)

## function for plotting the mean value
panel.mean <- function(x,y,...){
	x <- as.numeric(x)
	x.unique <- unique(x)
	for(X in x.unique) {
		Y <- y[x == X]
		if (!length(Y)) next
		mean.value <- list(y = mean(Y), x = X)
		do.call("lpoints", c(mean.value, pch = 20, col= "red"))
	}
}

# functions for the visualisation of different operatorprobabilities 
preprocess_op <- function(d, dname = "datasetname"){
	d <- as.data.frame(d)
    colnames(d) <- c("c0m50sp50", "c20m40sp40","c40m30sp30","c20m20sp60","c20m60sp20")
    x <- d*100
    x <- reshape(x, idvar="samp", times=names(x), timevar = "operatorprob",varying= list(names(x)), direction="long")
	names(x)[[2]] <- "value"
	x$ds <- dname
    return(x)
}

preprocess_op_comp <- function(ntrees, colum){	
	rt <- preprocess_op(d = rheart[,colum], dname = "Statlog heart")
	rt <- rbind(rt, preprocess_op(d = rcredit[,colum], dname = "Statlog German credit"))
	rt <- rbind(rt, preprocess_op(d = rspam[,colum], dname = "Spam"))
	rt <- rbind(rt, preprocess_op(d = rchessboard44_5[,colum], dname = "Chessboard 4x4 (5% noise)"))
	rt <- cbind(rt, rep( ntrees, dim(rt)[1]))
	colnames(rt)[5] <- "nIter"
	rt
}

r2 <- preprocess_op_comp("200 iterations", c(11:15))
r2 <- rbind(r2, preprocess_op_comp("500 iterations", c(16:20)))
r2 <- rbind(r2, preprocess_op_comp("10000 iterations", c(21:25)))


sort_op <- function(x){
	x$ds <- relevel(x$ds, "Statlog heart")
	x$ds <- relevel(x$ds, "Statlog German credit")
	x$ds <- relevel(x$ds, "Spam")
	x$ds <- relevel(x$ds, "Chessboard 4x4 (5% noise)")
	x
}

r2$operatorprob <- factor(r2$operatorprob)
r2$ds <- factor(r2$ds)
r2$nIter <- factor(r2$nIter)
r2 <- sort_op(r2)

b1 <- bwplot (value ~ factor(operatorprob)| nIter+ds , data= as.data.frame(r2), 
horizontal = FALSE,  
ylab= list("Accuracy (\\%)", cex=1.1),
xlab= list("Operator probability setting", cex=1.1),
pch= '|',
layout = c(3,4),
ylim=as.data.frame(matrix(c(
rep(c(55,100),3),
rep(c(85,94),3),
rep(c(60,80),3),
rep(c(58,90),3)
), nrow=2)),
scales= list(x= list(rot=60), y=list(relation="free"), alternating = F), 
panel=function(x,y,...) {
panel.bwplot(x,y,...)
panel.mean(x,y,...)
}
)

# functions for the visualisation of different population sizes 
preprocess_ntrees <- function(d, dname = "datasetname"){
	d <- as.data.frame(d)
    colnames(d) <- c("25 trees", "50 trees","100 trees","250 trees","500 trees")
    x <- d*100
    x <- reshape(x, idvar="samp", times=names(x), timevar = "ntrees",varying= list(names(x)), direction="long")
    names(x)[[2]] <- "value"
	x$ds <- dname
    return(x)
}

r <- preprocess_ntrees(d = rheart[,1:5], dname = "Statlog heart")
r <- rbind(r, preprocess_ntrees(d = rcredit[,1:5], dname = "Statlog German credit"))
r <- rbind(r, preprocess_ntrees(d = rspam[,1:5], dname = "Spam"))
r <- rbind(r, preprocess_ntrees(d = rchessboard44_5[,1:5], dname = "Chessboard 4x4 (5% noise)"))

sort_ntrees <- function(x){
	x$ds <- relevel(x$ds, "Chessboard 4x4 (5% noise)")
	x$ds <- relevel(x$ds, "Spam")
	x$ds <- relevel(x$ds, "Statlog German credit")
    x$ds <- relevel(x$ds, "Statlog heart")
	x$ntrees <- relevel(x$ntrees, "250 trees")
	x$ntrees <- relevel(x$ntrees, "100 trees")
	x$ntrees <- relevel(x$ntrees, "50 trees")
	x$ntrees <- relevel(x$ntrees, "25 trees")
	x
}

r$ntrees <- factor(r$ntrees) 
r$ds <- factor(r$ds)
r <- sort_ntrees(r) 
par.settings = list(cex=1.2)
b2 <- bwplot (value ~ ntrees | ds, data= as.data.frame(r), 
horizontal = FALSE,  
ylab= list("Accuracy (%)", cex=1.1),
pch= '|',
ylim=as.data.frame(matrix(c(
c(60,90),
c(65,80),
c(89,94),
c(60,95)
), nrow=2)),
scales= list(x= list(rot=60), y=list(relation="free"), alternating = FALSE), 
layout = c(4,1),
panel=function(x,y,...) {
        panel.bwplot(x,y,...)
        panel.mean(x,y,...)
       }
)


###################################################
### code chunk number 20: param-operator-plot
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
plot(b1)


###################################################
### code chunk number 21: param-ntrees-plot
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
plot(b2)


