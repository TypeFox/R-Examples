library(latticeExtra)
library(matrixTrellis)

data(apple, package="HH")
apple.aov.1 <- aov(yield ~ block + pre*treat, data=apple)
apple.aov.2 <- aov(yield ~ block + pre + treat, data=apple)
apple.aov.2b<- aov(yield ~ block + treat + pre, data=apple)

apple$yield.block.effect <- fitted(lm(yield ~ block, data=apple)) - mean(apple$yield)
apple$pre.block.effect   <- fitted(lm(pre   ~ block, data=apple)) - mean(apple$pre)
apple$yield.block        <- apple$yield - apple$yield.block.effect
apple$pre.block          <- apple$pre - apple$pre.block.effect


apple.ancova.3be<-  ancovaplot(yield.block ~ treat*pre.block, data=apple,
                               groups=block, pch=letters[1:4], cex=1.4,
                               col.line=trellis.par.get()$superpose.symbol$col,
                               col.by.groups=TRUE)
apple.ancova.3be

apple.ancova.3bh<-  ancovaplot(yield.block ~ treat*pre.block, data=apple,
                               col=c("red","green","purple","orange","blue","navyblue"),
                               pch=letters[1:6], cex=1.4,
                               col.line=trellis.par.get()$superpose.symbol$col,
                               col.by.groups=FALSE)
apple.ancova.3bh

apple.ancova.4 <-  ancovaplot(yield.block ~ pre.block + treat, data=apple)
apple.aov.4 <- aov(yield.block ~ pre.block + treat, data=apple)
apple.ancova.6 <-  ancovaplot(yield.block ~ treat, x=pre.block, data=apple)

apple$yield.block.pre <-
  apple$yield.block -
  predict.lm(apple.aov.4, type="terms", terms="pre.block")
apple.ancova.5 <-  ancovaplot(yield.block.pre ~ treat, x=pre.block, data=apple)
apple.ancova.7 <-  ancovaplot(yield.block ~ pre.block, groups=treat, data=apple)

apple.ancova.8b <-  ancovaplot(yield ~ pre * treat, data=apple,
                               groups=block, pch=letters[1:4], cex=1.4,
                               col.line=trellis.par.get()$superpose.symbol$col,
                               condition=apple$treat)
apple.ancova.8b

apple7 <-  do.call(c, list(a5 =apple.ancova.5,
                           a4 =apple.ancova.4,
                           a6 =apple.ancova.6,
                           a7 =apple.ancova.7,
                           a3be=apple.ancova.3be,
                           a8b=apple.ancova.8b,
                           layout=c(7, 6)))
apple7 <- update(apple7, scales=list(y=list(rot=0)))
apple7

apple7mt <- matrix.trellis(apple7, nrow=6, ncol=7, byrow=TRUE)
apple7mt$condlevels[[2]] <- c("5","4","6","7","3bd","8")
apple7mt$condlevels[[1]][7] <- "Superpose"
apple7mt

useOuterStrips(apple7mt)
combineLimits(apple7mt)
apple7uc <- useOuterStrips(combineLimits(apple7mt))
apple7uc
apple7ucr <- update(apple7uc, scales=list(
                                x=list(relation="same", at=7:11)))
apple7ucr

tmp <- apple7ucr$y.limits
dim(tmp) <- c(7,6)
tmp.range <- sapply(1:6, function(i) range(unlist(tmp[,i])))

for (i in 1:35) apple7ucr$y.limits[[i]] <- range(tmp.range[,-6])
for (i in 36:42) apple7ucr$y.limits[[i]] <- tmp.range[,6]
apple7ucr <- update(apple7ucr,
                    xlab="pre (adjusted as indicated)",
                    ylab="yield (adjusted as indicated)",
                    scales=list(
                      y=list(at=rep(list(c(200,250,300,350),
                               NULL,NULL,NULL,NULL,NULL,NULL), 6)),
                      x=list(alternating=1)),
                    xlab.top="Treatment")
apple7ucr

apple7ucrs <-
  resizePanels(update(apple7ucr,
                      between=list(
                        x=c(.2,.2,.2,.2,.2,1),
                        y=c(.5,.2,.2,.2,1))),
               h=c(rep(diff(apple7ucr$y.limits[[1]]), 5),
                 diff(apple7ucr$y.limits[[36]])))
apple7ucrs

apple7ucrs2 <- apple7ucrs
apple7ucrs2$condlevels[[2]] <-
  c("y|bx~t","y|b~x|b+t","y|b~t","y|b~x|b","y|b~x|b*t","y~x*t")
apple7ucrs2
