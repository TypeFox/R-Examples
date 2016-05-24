### R code from vignette source 'ShortIntro.Rnw'

###################################################
### code chunk number 1: Init_hidden
###################################################
library(boilerpipeR)
data(content)
options(width = 60)


###################################################
### code chunk number 2: ShortIntro.Rnw:56-58 (eval = FALSE)
###################################################
## library(boilerpipeR)
## library(RCurl)


###################################################
### code chunk number 3: ShortIntro.Rnw:61-63 (eval = FALSE)
###################################################
## url <- "http://quantivity.wordpress.com/2012/11/09/multi-asset-market-regimes/"
## content <- getURL(url)


###################################################
### code chunk number 4: ShortIntro.Rnw:71-72
###################################################
cat(substr(content, 1, 80))


###################################################
### code chunk number 5: ShortIntro.Rnw:94-96
###################################################
extract <- DefaultExtractor(content)
cat(substr(extract, 1, 120))


###################################################
### code chunk number 6: ShortIntro.Rnw:120-127
###################################################
articleextract <- ArticleExtractor(content)
articlesentencesextract <- ArticleSentencesExtractor(content)
canolaextract <- CanolaExtractor(content)
defaultextract <- DefaultExtractor(content)
keepeverythingextract <- KeepEverythingExtractor(content)
largestcontentextract <- LargestContentExtractor(content)
numwordsrulesextract <- NumWordsRulesExtractor(content)


