### R code from vignette source 'RcppDE.Rnw'

###################################################
### code chunk number 1: prelim
###################################################
options(width=50)
library(lattice)
library(RcppDE)
RcppDE.version <- packageDescription("RcppDE")$Version
RcppDE.date <- packageDescription("RcppDE")$Date
now.date <- strftime(Sys.Date(), "%B %d, %Y")
# create figures/ if not present
if ( ! (file.exists("figures") && file.info("figures")[["isdir"]]) ) dir.create("figures")


###################################################
### code chunk number 2: smallRes
###################################################
## # small benchmark at SVN 2419M
## # At 2010-11-08 06:42:29.018531
smallLines <- "
            DEoptim   RcppDE ratioRcppToBasic pctGainOfRcpp netSpeedUp
Rastrigin5  0.10912 0.099875          0.91523        8.4765     1.0926
Rastrigin10 0.23738 0.214875          0.90521        9.4787     1.1047
Rastrigin20 0.55587 0.501500          0.90218        9.7819     1.1084
Wild5       0.18288 0.171875          0.93985        6.0150     1.0640
Wild10      0.40912 0.391125          0.95600        4.3996     1.0460
Wild20      1.04513 0.987375          0.94474        5.5257     1.0585
Genrose5    0.18913 0.179250          0.94779        5.2214     1.0551
Genrose10   0.39538 0.374625          0.94752        5.2482     1.0554
Genrose20   0.90050 0.848375          0.94212        5.7885     1.0614
"
## MEANS       0.44717 0.418764          0.93648        6.3517     1.0678
## # Done 2010-11-08 06:43:50.88171

con <- textConnection(smallLines)
smallData <- read.table(con, header=TRUE, sep="")
close(con)

sb <- trellis.par.get("strip.background")
sb[["col"]][1:2] <- c("gray80","gray90")
trellis.par.set("strip.background", sb)

# dput(brewer.pal(7, "Set1"))
.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
           "#FF7F00", "#FFF33", "#A65628")[-6]

ss <- trellis.par.get("superpose.symbol")
ss[["col"]][1:6] <- .cols
ss[["cex"]] <- rep(1.0, 7)
ss[["pch"]] <- rep(19, 7)
ss[["alpha"]] <- rep(0.75, 7)
trellis.par.set("superpose.symbol", ss)

smallWide <- data.frame(timeInSeconds=c(smallData[,1,drop=TRUE], smallData[,2,drop=TRUE]),
                        pkg=rep(c("DEoptim", "RcppDE"), each=9),
                        fun=rep(rep(c("Rastrigin", "Wild", "Genrose"), each=3), 2),
                        n=c(5,10,20), 6)
smallWide$fun <- factor(smallWide$fun, levels=c("Rastrigin", "Genrose", "Wild"))
print(dotplot(as.factor(n) ~ timeInSeconds | fun, group=pkg, data=smallWide, layout=c(1,3),
              xlab="Time in seconds for 5, 10 and 20 parameter problems, using logarithmic axis", ylab="",
              scales=list(x=list(log=TRUE,at=c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0), labels=c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0))),
              key=simpleKey(text=c("DEoptim","RcppDE"), space="top")))


###################################################
### code chunk number 3: largeRes
###################################################
## # big benchmark at SVN 2419M
## # At 2010-11-08 06:43:51.422299
largeLines <- "
             DEoptim RcppDE ratioRcppToBasic pctGainOfRcpp netSpeedUp
Rastrigin50    1.770  1.575          0.88983       11.0169     1.1238
Rastrigin100   4.794  4.258          0.88819       11.1806     1.1259
Rastrigin200  14.840 12.472          0.84043       15.9569     1.1899
Wild50         3.692  3.558          0.96371        3.6295     1.0377
Wild100       11.127 10.646          0.95677        4.3228     1.0452
Wild200       38.026 35.755          0.94028        5.9722     1.0635
Genrose50      2.587  2.414          0.93313        6.6873     1.0717
Genrose100     6.252  5.739          0.91795        8.2054     1.0894
Genrose200    17.058 15.147          0.88797       11.2030     1.1262
"
## MEANS         11.127 10.174          0.91431        8.5695     1.0937
## # Done 2010-11-08 06:47:03.810348

con <- textConnection(largeLines)
largeData <- read.table(con, header=TRUE, sep="")
close(con)

largeWide <- data.frame(timeInSeconds=c(largeData[,1,drop=TRUE], largeData[,2,drop=TRUE]),
                        pkg=rep(c("DEoptim", "RcppDE"), each=9),
                        fun=rep(rep(c("Rastrigin", "Wild", "Genrose"), each=3), 2),
                        n=c(50,100,200), 6)
largeWide$fun <- factor(largeWide$fun, levels=c("Rastrigin", "Genrose", "Wild"))
print(dotplot(as.factor(n) ~ timeInSeconds | fun, group=pkg, data=largeWide, layout=c(1,3),
              xlab="Time in seconds for 50, 100 and 200 parameter problems, using logarithmic axis", ylab="",
              scales=list(x=list(log=TRUE,at=c(1, 2, 5, 10, 20, 30), labels=c(1, 2, 5, 10, 20, 30))),
              key=simpleKey(text=c("DEoptim","RcppDE"), space="top")))


###################################################
### code chunk number 4: compiledRes
###################################################
## # compiled benchmark at SVN 2419:2421M
## # At 2010-11-08 06:48:42.56918
compiledLines <- "
             DEoptim  RcppDE ratioRcppToBasic pctGainOfRcpp netSpeedUp
Rastrigin50    1.781  0.6090          0.34194        65.806     2.9245
Rastrigin100   4.807  2.0940          0.43561        56.439     2.2956
Rastrigin200  14.572  7.5000          0.51469        48.531     1.9429
Wild50         3.748  0.9500          0.25347        74.653     3.9453
Wild100       11.268  3.3160          0.29428        70.572     3.3981
Wild200       37.225 12.4120          0.33343        66.657     2.9991
Genrose50      2.667  0.2710          0.10161        89.839     9.8413
Genrose100     6.498  0.7190          0.11065        88.935     9.0376
Genrose200    17.471  1.9830          0.11350        88.650     8.8104
"
## MEANS         11.115  3.3171          0.29843        70.157     3.3509
## # Done 2010-11-08 06:50:53.195003

con <- textConnection(compiledLines)
compiledData <- read.table(con, header=TRUE, sep="")
close(con)

compiledWide <- data.frame(timeInSeconds=c(compiledData[,1,drop=TRUE], compiledData[,2,drop=TRUE]),
                        pkg=rep(c("DEoptim", "RcppDE"), each=9),
                        fun=rep(rep(c("Rastrigin", "Wild", "Genrose"), each=3), 2),
                        n=c(50,100,200), 6)
compiledWide$fun <- factor(compiledWide$fun, levels=c("Rastrigin", "Genrose", "Wild"))
print(dotplot(as.factor(n) ~ timeInSeconds | fun, group=pkg, data=compiledWide, layout=c(1,3),
              xlab="Time in sec. for 50, 100 and 200 parameter problems, compiled objective function, logarithmic axis", ylab="",
              scales=list(x=list(log=TRUE,at=c(0.5, 1, 2, 5, 10, 20, 30), labels=c(0.5, 1, 2, 5, 10, 20, 30))),
              key=simpleKey(text=c("DEoptim","RcppDE"), space="top")))


