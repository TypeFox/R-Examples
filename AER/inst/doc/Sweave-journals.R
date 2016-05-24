### R code from vignette source 'Sweave-journals.Rnw'

###################################################
### code chunk number 1: Sweave-journals.Rnw:8-11
###################################################
data("Journals", package = "AER")
journals_lm <- lm(log(subs) ~ log(price/citations), data = Journals)
journals_lm


###################################################
### code chunk number 2: Sweave-journals.Rnw:17-19
###################################################
plot(log(subs) ~ log(price/citations), data = Journals)
abline(journals_lm)


