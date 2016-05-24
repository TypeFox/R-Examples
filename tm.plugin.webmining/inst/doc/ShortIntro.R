### R code from vignette source 'ShortIntro.Rnw'

###################################################
### code chunk number 1: Init_hidden
###################################################
library(tm)
library(tm.plugin.webmining)
data(yahoonews)
options(width = 60)


###################################################
### code chunk number 2: ShortIntro.Rnw:55-57 (eval = FALSE)
###################################################
## library(tm)
## library(tm.plugin.webmining)


###################################################
### code chunk number 3: ShortIntro.Rnw:68-69 (eval = FALSE)
###################################################
## yahoonews <- WebCorpus(YahooNewsSource("Microsoft"))


###################################################
### code chunk number 4: ShortIntro.Rnw:78-79
###################################################
class(yahoonews)


###################################################
### code chunk number 5: ShortIntro.Rnw:86-87
###################################################
yahoonews


###################################################
### code chunk number 6: ShortIntro.Rnw:108-115
###################################################
# Little hack to restrict output width
meta(yahoonews[[1]], "description") <- 
		paste(substring(meta(yahoonews[[1]], "description"), 1, 70), "...", sep = "")
meta(yahoonews[[1]], "id") <- 
		paste(substring(meta(yahoonews[[1]], "id"), 1, 70), "...", sep = "")
meta(yahoonews[[1]], "origin") <- 
		paste(substring(meta(yahoonews[[1]], "origin"), 1, 70), "...", sep = "")


###################################################
### code chunk number 7: ShortIntro.Rnw:117-118
###################################################
meta(yahoonews[[1]])


###################################################
### code chunk number 8: ShortIntro.Rnw:126-129
###################################################
# Little hack to restrict output length
content(yahoonews[[1]]) <- 
		paste(substring(yahoonews[[1]], 1, 100), "...", sep = "")


###################################################
### code chunk number 9: ShortIntro.Rnw:131-132
###################################################
yahoonews[[1]]


###################################################
### code chunk number 10: ShortIntro.Rnw:139-140 (eval = FALSE)
###################################################
## inspect(yahoonews)


###################################################
### code chunk number 11: ShortIntro.Rnw:159-166 (eval = FALSE)
###################################################
## googlefinance <- WebCorpus(GoogleFinanceSource("NASDAQ:MSFT"))
## googlenews <- WebCorpus(GoogleNewsSource("Microsoft"))
## nytimes <- WebCorpus(NYTimesSource("Microsoft", appid = nytimes_appid))
## reutersnews <- WebCorpus(ReutersNewsSource("businessNews"))
## yahoofinance <- WebCorpus(YahooFinanceSource("MSFT"))
## yahooinplay <- WebCorpus(YahooInplaySource())
## yahoonews <- WebCorpus(YahooNewsSource("Microsoft"))


###################################################
### code chunk number 12: ShortIntro.Rnw:180-181 (eval = FALSE)
###################################################
## yahoonews <- corpus.update(yahoonews)


