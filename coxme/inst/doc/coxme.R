### R code from vignette source 'coxme.Rnw'

###################################################
### code chunk number 1: coxme.Rnw:110-111 (eval = FALSE)
###################################################
## fit1 <- coxme(Surv(endage, cancer) ~ parity + (1| famid))


###################################################
### code chunk number 2: coxme.Rnw:142-147
###################################################
library(coxme)
stem(table(eortc$center))

efit1 <- coxph(Surv(y, uncens) ~ trt, eortc)
efit2 <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc)


###################################################
### code chunk number 3: coxme.Rnw:153-154
###################################################
print(efit2)


###################################################
### code chunk number 4: coxme.Rnw:241-242
###################################################
stem(exp(ranef(efit2)[[1]]))


###################################################
### code chunk number 5: coxme.Rnw:247-249
###################################################
efit3 <- coxme(Surv(y, uncens) ~ trt + (1 | center/trt), eortc)
efit3


###################################################
### code chunk number 6: coxme.Rnw:292-309
###################################################
library(coxme)
library(kinship2)
options(show.signif.stars=FALSE)
makefig <- function(file) {
    pdf(paste(file, "pdf", sep='.'), width=7, height=5)
    par(mar=c(5.1, 4.1, .1, .1))
}

names(minnbreast)
with(minnbreast, table(sex, cancer, exclude=NULL))

mped <- with(minnbreast, pedigree(id, fatherid, motherid, sex,
                                  affected=cancer, famid=famid,
                                  status=proband))
makefig("cfig1")
plot(mped["8"])  #figure 1
dev.off()


###################################################
### code chunk number 7: coxme.Rnw:334-342
###################################################
minnfemale <- minnbreast[minnbreast$sex == 'F' & !is.na(minnbreast$sex),]
fit0 <- coxph(Surv(endage, cancer) ~ I(parity>0), minnfemale,
              subset=(proband==0))
summary(fit0)

fit1 <- coxme(Surv(endage, cancer) ~ I(parity>0) + (1|famid),
              minnfemale, subset=(proband==0))
print(fit1)


###################################################
### code chunk number 8: coxme.Rnw:380-399
###################################################
ncancer <- with(minnfemale, tapply(cancer, famid, sum, na.rm=T))
pyears <-  with(minnfemale, tapply(endage -18, famid, sum, na.rm=T))
count  <-  with(minnfemale, tapply(cancer, famid, 
                                   function(x) sum(!is.na(x))))
indx <- match(names(ranef(fit1)[[1]]), names(ncancer))                

makefig("cfig2")
plot(ncancer[indx], exp(ranef(fit1)[[1]]), log='y',
     xlab="Number of cancers per family",
     ylab="Estimated familial risk")
abline(h=1, lty=2)
text(c(8.1, 1.6), c(.85, 1.2), c("165", "72"))
dev.off()

indx <- match(c(72,165), names(ncancer))
temp <- cbind(ncancer, count, pyears, 100*ncancer/pyears)[indx,]
dimnames(temp) <- list(c(72, 165), 
                       c("Cancers", "N", "Years of FU", "Rate"))
print(round(temp,2))


###################################################
### code chunk number 9: coxme.Rnw:423-436
###################################################
estvar <- seq(.2, .6, length=15)^2  #range of std values
loglik <- double(15)
for (i in 1:15) {
    tfit <- coxme(Surv(endage, cancer) ~ I(parity>0) + (1|famid),
                  data=minnfemale, subset=(proband==0),
                  vfixed=estvar[i])
    loglik[i] <- 2*diff(tfit$loglik)[1]
}
makefig("cfig3")
plot(sqrt(estvar), loglik, 
     xlab="Std of the random effect", ylab="2 * loglik")
abline(h=2*diff(fit1$loglik)[1] - qchisq(.95, 1), lty=2)
dev.off()


###################################################
### code chunk number 10: coxme.Rnw:447-450
###################################################
temp <- 2*diff(fit1$loglik)[1] - loglik
approx(temp[1:8], sqrt(estvar[1:8]), 3.84)$y
approx(temp[9:15], sqrt(estvar[9:15]), 3.84)$y


###################################################
### code chunk number 11: coxme.Rnw:488-493
###################################################
kmat <- kinship(mped)
fit2 <- coxme(Surv(endage, cancer) ~ I(parity>0) + (1|id),
              data=minnfemale, varlist=coxmeMlist(2*kmat, rescale=F),
              subset=(proband==0))
print(fit2)


