### R code from vignette source 'Using_Hot_Deck_Data.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(useFancyQuotes=FALSE, width=100)


###################################################
### code chunk number 2: Using_Hot_Deck_Data.Rnw:36-39
###################################################
library(hot.deck)
data(isq99)
out <- hot.deck(isq99, sdCutoff=3, IDvars = c("IDORIGIN", "YEAR"))


###################################################
### code chunk number 3: numdonors
###################################################
numdonors <- sapply(out$donors, length)
numdonors <- sapply(out$donors, length)
numdonors <- ifelse(numdonors > 5, 6, numdonors)
numdonors <- factor(numdonors, levels=1:6, labels=c(1:5, ">5"))
table(numdonors)


###################################################
### code chunk number 4: tscslag
###################################################
tscslag <- function(dat, x, id, time){
	obs <- apply(dat[, c(id, time)], 1, paste, collapse=".")
	tm1 <- dat[[time]] - 1
	lagobs <- apply(cbind(dat[[id]], tm1), 1, paste, collapse=".")
	lagx <- dat[match(lagobs, obs), x]
}
for(i in 1:length(out$data)){
    out$data[[i]]$lagAI <- tscslag(out$data[[i]], "AI", "IDORIGIN", "YEAR")
    out$data[[i]]$lagPCGNP <- tscslag(out$data[[i]], "PCGNP", "IDORIGIN", "YEAR")
    out$data[[i]]$lagLPOP <- tscslag(out$data[[i]], "LPOP", "IDORIGIN", "YEAR")
}


###################################################
### code chunk number 5: pcgchange
###################################################
for(i in 1:length(out$data)){
    out$data[[i]]$pctchgPCGNP <- with(out$data[[i]], c(PCGNP-lagPCGNP)/lagPCGNP)
    out$data[[i]]$pctchgLPOP <- with(out$data[[i]], c(LPOP-lagLPOP)/lagLPOP)
}


###################################################
### code chunk number 6: hd2amelia
###################################################
out <- hd2amelia(out)


###################################################
### code chunk number 7: zel
###################################################
library(Zelig)


###################################################
### code chunk number 8: zelig
###################################################
library(Zelig)
z <- zelig(AI ~ lagAI + pctchgPCGNP + PCGNP + pctchgLPOP + LPOP + MIL2 + LEFT + 
    BRIT + POLRT + CWARCOW + IWARCOW2, data=out, model="normal", cite=FALSE)
summary(z)


###################################################
### code chunk number 9: mods
###################################################
# initialize list
results <- list()
# loop over imputed datasets
for(i in 1:length(out$imputations)){
    results[[i]] <- lm(AI ~ lagAI + pctchgPCGNP + PCGNP + pctchgLPOP + LPOP + MIL2 + LEFT + 
    BRIT + POLRT + CWARCOW + IWARCOW2, data=out$imputations[[i]])
}
library(mitools)
summary(MIcombine(results))


###################################################
### code chunk number 10: conv
###################################################
library(miceadds)
out.mids <- datalist2mids(out$imputations)
s <- summary(pool(lm.mids(AI ~ lagAI + pctchgPCGNP + PCGNP + pctchgLPOP + LPOP + MIL2 + LEFT + 
BRIT + POLRT + CWARCOW + IWARCOW2, data=out.mids)))
round(s, 4)


