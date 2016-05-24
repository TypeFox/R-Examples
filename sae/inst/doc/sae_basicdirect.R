### R code from vignette source 'sae_basicdirect.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: sae_basicdirect.Rnw:139-143
###################################################
library("sae")
data("incomedata")
data("sizeprov")
data("sizeprovedu")


###################################################
### code chunk number 2: sae_basicdirect.Rnw:152-154
###################################################
z <- 6557.143
poor <- as.integer(incomedata$income < z)


###################################################
### code chunk number 3: sae_basicdirect.Rnw:160-163
###################################################
Popn <- sizeprov[, c("provlab", "Nd")]
DIR <- direct(y = poor, dom = incomedata$provlab, 
              sweight = incomedata$weight, domsize = Popn)


###################################################
### code chunk number 4: sae_basicdirect.Rnw:176-181
###################################################
Popn.educ <- sizeprovedu[,-2]
colnames(Popn.educ) <- c("provlab","0", "1", "2", "3")
PSYN.educ <- pssynt(y = poor, sweight = incomedata$weight, 
                    ps = incomedata$educ, 
                    domsizebyps = Popn.educ)


###################################################
### code chunk number 5: sae_basicdirect.Rnw:190-193
###################################################
SSD <- ssd(dom = provlab, sweight = weight, domsize = Popn, 
           direct = DIR[, c("Domain", "Direct")], 
           synthetic = PSYN.educ, data = incomedata)


###################################################
### code chunk number 6: sae_basicdirect.Rnw:198-204
###################################################
results <- data.frame(Province = DIR$Domain, 
                      SampleSize = DIR$SampSize, 
                      DIR = DIR$Direct * 100, 
                      PSYN.educ = PSYN.educ$PsSynthetic * 100, 
                      SSD = SSD$ssd * 100)
print(results, row.names = FALSE)


###################################################
### code chunk number 7: sae_basicdirect.Rnw:221-233
###################################################
# Sorted results by decreasing sample size
results <- results[order(results$SampleSize, 
                         decreasing = TRUE), ]
plot(results$DIR, type = "n", 
     xlab = "area (sorted by decreasing sample size)", 
     ylab = "Estimate", cex.axis = 1.5, cex.lab = 1.5)
points(results$DIR, type = "b", col = 3, lwd = 2, pch = 1)
points(results$PSYN.educ, type= "b", col = 5, lwd = 2, pch = 2)
points(results$SSD, type = "b", col = 2, lwd = 2, pch = 5)
legend("bottom", legend = c("Direct", "Post-strat educ", "SSD"), 
       ncol = 1, col = c(3, 5, 2), lwd = rep(2, 3), 
       pch = c(1, 2, 5), cex = 1.3)


###################################################
### code chunk number 8: sae_basicdirect.Rnw:240-241
###################################################
DIR[c("42","5","34","44","40"), -4]


