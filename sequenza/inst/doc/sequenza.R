### R code from vignette source 'sequenza.Rnw'

###################################################
### code chunk number 1: sequenza.Rnw:60-61
###################################################
  options(width = 60, strict.width = "wrap")


###################################################
### code chunk number 2: instLib (eval = FALSE)
###################################################
## install.packages("sequenza")


###################################################
### code chunk number 3: findexec (eval = FALSE)
###################################################
## system.file("exec", "sequenza-utils.py", package="sequenza")


###################################################
### code chunk number 4: varscan (eval = FALSE)
###################################################
## snp <- read.table("varscan.snp", header = TRUE, sep = "\t")
## cnv <- read.table("varscan.copynumber", header = TRUE, sep = "\t")
## seqz.data <- VarScan2seqz(varscan.somatic = snp, varscan.copynumber = cnv)
## 
## write.table(seqz.data, "my.sample.seqz", col.names = TRUE, row.names = FALSE, sep = "\t")


###################################################
### code chunk number 5: loadLib
###################################################
library("sequenza")


###################################################
### code chunk number 6: setFile (eval = FALSE)
###################################################
## data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")
## data.file


###################################################
### code chunk number 7: setFile2
###################################################
data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")


###################################################
### code chunk number 8: readAfreqChr (eval = FALSE)
###################################################
## seqz.data <- read.seqz(data.file, chr.name = "1")


###################################################
### code chunk number 9: readAfreq
###################################################
seqz.data <- read.seqz(data.file)


###################################################
### code chunk number 10: sequenza.Rnw:202-203
###################################################
str(seqz.data, vec.len = 2)


###################################################
### code chunk number 11: depthRGCnormall
###################################################
gc.stats <- gc.sample.stats(data.file)


###################################################
### code chunk number 12: gcstr
###################################################
str(gc.stats)


###################################################
### code chunk number 13: depthRGCnorm (eval = FALSE)
###################################################
## gc.stats <- gc.norm(x = seqz.data$depth.ratio,
##                     gc = seqz.data$GC.percent)


###################################################
### code chunk number 14: useGCmedians
###################################################
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)

seqz.data$adjusted.ratio <- seqz.data$depth.ratio / 
                           gc.vect[as.character(seqz.data$GC.percent)]                       


###################################################
### code chunk number 15: GCfig
###################################################
par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
matplot(gc.stats$gc.values, gc.stats$raw,
        type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2),
        xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio,
      breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2),
      xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')


###################################################
### code chunk number 16: sequenzaExtract
###################################################
test <- sequenza.extract(data.file)
names(test)


###################################################
### code chunk number 17: 3panelsPlot
###################################################

chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]], 
                ratio.windows = test$ratio[[1]], min.N.ratio = 1,
                segments = test$segments[[1]], main = test$chromosomes[1])


###################################################
### code chunk number 18: sequenzaFit (eval = FALSE)
###################################################
## CP.example <- sequenza.fit(test)


###################################################
### code chunk number 19: loadCP
###################################################
data(CP.example)


###################################################
### code chunk number 20: sequenzaRes (eval = FALSE)
###################################################
## sequenza.results(sequenza.extract = test, cp.table = CP.example,
##                  sample.id = "Test", out.dir="TEST")


###################################################
### code chunk number 21: ConfIntCP
###################################################
cint <- get.ci(CP.example)


###################################################
### code chunk number 22: CPplot
###################################################
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE, likThresh = c(0.95))


###################################################
### code chunk number 23: CPplotCI
###################################################
par(mfrow = c(2,2))
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE)
plot(cint$values.cellularity, ylab = "Cellularity",
     xlab = "posterior probability", type = "n")
select <- cint$confint.cellularity[1] <= cint$values.cellularity[,2] &
          cint$values.cellularity[,2] <= cint$confint.cellularity[2]
polygon(y = c(cint$confint.cellularity[1], cint$values.cellularity[select, 2], cint$confint.cellularity[2]), 
        x = c(0, cint$values.cellularity[select, 1], 0), col='red', border=NA)
lines(cint$values.cellularity)
abline(h = cint$max.cellularity, lty = 2, lwd = 0.5)  

plot(cint$values.ploidy, xlab = "Ploidy",
     ylab = "posterior probability", type = "n")
select <- cint$confint.ploidy[1] <= cint$values.ploidy[,1] &
          cint$values.ploidy[,1] <= cint$confint.ploidy[2]
polygon(x = c(cint$confint.ploidy[1], cint$values.ploidy[select, 1], cint$confint.ploidy[2]), 
        y = c(0, cint$values.ploidy[select, 2], 0), col='red', border=NA)
lines(cint$values.ploidy)
abline(v = cint$max.ploidy, lty = 2, lwd = 0.5)



###################################################
### code chunk number 24: seParam
###################################################
cellularity <- cint$max.cellularity
cellularity
ploidy <- cint$max.ploidy
ploidy


###################################################
### code chunk number 25: avgDepth
###################################################
avg.depth.ratio <- mean(test$gc$adj[, 2])
avg.depth.ratio


###################################################
### code chunk number 26: mmufBayes
###################################################
mut.tab     <- na.exclude(do.call(rbind, test$mutations))
mut.alleles <- mufreq.bayes(mufreq = mut.tab$F,
                            depth.ratio = mut.tab$adjusted.ratio,
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio)

head(mut.alleles)
head(cbind(mut.tab[,c("chromosome","position","F","adjusted.ratio", "mutation")],
           mut.alleles))


###################################################
### code chunk number 27: bafBayes
###################################################
seg.tab     <- na.exclude(do.call(rbind, test$segments))
cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio,
                        cellularity = cellularity, ploidy = ploidy,
                        avg.depth.ratio = avg.depth.ratio)
head(cn.alleles)
seg.tab <- cbind(seg.tab, cn.alleles)
head(seg.tab)


###################################################
### code chunk number 28: chrViewWithCP
###################################################
chromosome.view(mut.tab = test$mutations[[3]], baf.windows = test$BAF[[3]], 
                ratio.windows = test$ratio[[3]],  min.N.ratio = 1,
                segments = seg.tab[seg.tab$chromosome == test$chromosomes[3],],
                main = test$chromosomes[3],
                cellularity = cellularity, ploidy = ploidy,
                avg.depth.ratio = avg.depth.ratio)


###################################################
### code chunk number 29: genomeViewCNt
###################################################
genome.view(seg.cn = seg.tab, info.type = "CNt")
legend("bottomright", bty="n", c("Tumor copy number"),col = c("red"), 
       inset = c(0, -0.4), pch=15, xpd = TRUE)


###################################################
### code chunk number 30: genomeViewAB
###################################################
genome.view(seg.cn = seg.tab, info.type = "AB")
legend("bottomright", bty = "n", c("A-allele","B-allele"), col= c("red", "blue"), 
       inset = c(0, -0.45), pch = 15, xpd = TRUE)


