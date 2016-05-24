### R code from vignette source 'tileHMM.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: tileHMM.Rnw:38-41
###################################################
library(tileHMM)
data(simChIP)



###################################################
### code chunk number 2: tileHMM.Rnw:60-62
###################################################
h3k27.data <- as.matrix(simChIP[ , 7:14])



###################################################
### code chunk number 3: tileHMM.Rnw:67-71
###################################################

h3k27.data <- log(h3k27.data, 2)
h3k27.norm <- preprocessCore::normalize.quantiles(h3k27.data)



###################################################
### code chunk number 4: tileHMM.Rnw:78-88
###################################################
par(mfrow = c(2, 2))
smoothScatter(h3k27.norm[ , 1], h3k27.norm[ , 5], xlab = "H3 sample 1", 
    ylab = "H3K27me3 sample 1")
smoothScatter(h3k27.norm[ , 2], h3k27.norm[ , 6], xlab = "H3 sample 2", 
    ylab = "H3K27me3 sample 2")
smoothScatter(h3k27.norm[ , 3], h3k27.norm[ , 7], xlab = "H3 sample 3", 
    ylab = "H3K27me3 sample 3")
smoothScatter(h3k27.norm[ , 4], h3k27.norm[ , 8], xlab = "H3 sample 4", 
    ylab = "H3K27me3 sample 4")



###################################################
### code chunk number 5: tileHMM.Rnw:103-105
###################################################
h3k27.stat <- shrinkt.st(h3k27.norm, c(rep(2, 4), rep(1, 4)), verbose = 0)



###################################################
### code chunk number 6: tileHMM.Rnw:113-124
###################################################
def.par <- par(no.readonly = TRUE)

h <- hist(h3k27.stat, breaks = 100, plot=FALSE)
layout(matrix(c(1,2), 1, 2), width=c(3,1))
par(mar=c(5,4,2,1))
smoothScatter(simChIP[,2], h3k27.stat, xlab = 'genomic position', 
    ylab='probe statistic',nrpoints=500) 
par(mar=c(5,0,2,1))
barplot(h$counts, horiz=TRUE, axes=FALSE)

par(def.par)


###################################################
### code chunk number 7: tileHMM.Rnw:146-149
###################################################
hmm.init <- hmm.setup(h3k27.stat, 
    state = c("non-enriched","enriched"), pos.state = 2)



###################################################
### code chunk number 8: tileHMM.Rnw:155-158
###################################################
show(hmm.init)
plot(hmm.init)



###################################################
### code chunk number 9: tileHMM.Rnw:162-163
###################################################
plot(hmm.init)


###################################################
### code chunk number 10: tileHMM.Rnw:180-184
###################################################
max.gap <- 999
h3k27.gap <- diff(simChIP[['position']]) > max.gap
gap.idx <- which(h3k27.gap)



###################################################
### code chunk number 11: tileHMM.Rnw:191-196
###################################################
start <- c(1, gap.idx + 1)
end <- c(gap.idx, length(h3k27.stat))
h3k27.lst <- mapply(function(s, e, data) data[s:e], start, end, 
    MoreArgs = list(h3k27.stat)) 



###################################################
### code chunk number 12: tileHMM.Rnw:204-206
###################################################
hmm.opt <- viterbiEM(hmm.init, h3k27.lst, df = 9)



###################################################
### code chunk number 13: tileHMM.Rnw:245-248
###################################################
print(hmm.opt)
plot(hmm.opt)



###################################################
### code chunk number 14: tileHMM.Rnw:252-253
###################################################
plot(hmm.opt)


###################################################
### code chunk number 15: tileHMM.Rnw:264-268
###################################################
post <- lapply(h3k27.lst, posterior, hmm.opt)
state.seq <- lapply(post, apply, 2, which.max)
state.seq <- states(hmm.opt)[c(state.seq, recursive=TRUE)]



###################################################
### code chunk number 16: densOverlay (eval = FALSE)
###################################################
## ratio <- sum(state.seq == 'enriched') / length(state.seq)
## 
## hist(h3k27.stat, breaks=100, probability=TRUE, main='', 
##     xlab = 'probe statistic')
##     
## plot(hmm.opt@emission$enriched, new.plot = FALSE, 
##     lty = 1, lwd = 2, col = 4, weight = ratio)
## 
## plot(hmm.opt@emission$'non-enriched', new.plot = FALSE, 
##     lty = 2, lwd = 2, col = 2, weight = 1 - ratio)
## 
## legend('topright', legend = c('enriched', 'non-enriched'), 
##     lty = 1:2, lwd = 2, col = c(4,2))
## 


###################################################
### code chunk number 17: tileHMM.Rnw:298-299
###################################################
ratio <- sum(state.seq == 'enriched') / length(state.seq)

hist(h3k27.stat, breaks=100, probability=TRUE, main='', 
    xlab = 'probe statistic')
    
plot(hmm.opt@emission$enriched, new.plot = FALSE, 
    lty = 1, lwd = 2, col = 4, weight = ratio)

plot(hmm.opt@emission$'non-enriched', new.plot = FALSE, 
    lty = 2, lwd = 2, col = 2, weight = 1 - ratio)

legend('topright', legend = c('enriched', 'non-enriched'), 
    lty = 1:2, lwd = 2, col = c(4,2))



###################################################
### code chunk number 18: tileHMM.Rnw:312-314
###################################################
regions.idx <- region.position(state.seq, region = 'enriched')



###################################################
### code chunk number 19: tileHMM.Rnw:324-327
###################################################
regions.pos <- matrix(simChIP[regions.idx, 2], 
    nrow = 2, ncol = dim(regions.idx)[2])



###################################################
### code chunk number 20: tileHMM.Rnw:348-353
###################################################
post.enriched <- lapply(post,"[",2,)
post.enriched <- exp(c(post.enriched, recursive = TRUE))
region.score <- apply(regions.idx, 2, 
    function(reg, post) mean(post[reg[1]:reg[2]]), post.enriched)



###################################################
### code chunk number 21: tileHMM.Rnw:357-359
###################################################
region.len <- apply(regions.pos, 2, diff)



###################################################
### code chunk number 22: tileHMM.Rnw:369-375
###################################################
plot(region.len, region.score, xlab = 'region length',
    ylab = 'score', pch = '')
plot.idx <- region.len < 400 & region.score < 0.8
points(region.len[plot.idx], region.score[plot.idx], col = 2, pch = 4)
points(region.len[!plot.idx], region.score[!plot.idx], col = 4)



###################################################
### code chunk number 23: tileHMM.Rnw:383-386
###################################################
regions.clean <- remove.short(regions.idx, post.enriched, 
    simChIP[ , 1:2], min.length = 400, min.score = 0.8)



###################################################
### code chunk number 24: tileHMM.Rnw:399-400
###################################################
gff <- reg2gff(regions.clean, post.enriched, simChIP[ , 1:2])


###################################################
### code chunk number 25: tileHMM.Rnw:402-404 (eval = FALSE)
###################################################
## gff <- reg2gff(regions.clean, post.enriched, simChIP[ , 1:2],
##     file = 'simRegions.gff')


