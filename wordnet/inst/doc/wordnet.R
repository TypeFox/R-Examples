### R code from vignette source 'wordnet.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: wordnet.Rnw:15-16
###################################################
options(width = 75)


###################################################
### code chunk number 2: wordnet.Rnw:46-47
###################################################
library("wordnet")


###################################################
### code chunk number 3: wordnet.Rnw:66-67
###################################################
getFilterTypes()


###################################################
### code chunk number 4: wordnet.Rnw:77-80 (eval = FALSE)
###################################################
## filter <- getTermFilter("StartsWithFilter", "car", TRUE)
## terms <- getIndexTerms("NOUN", 5, filter)
## sapply(terms, getLemma)


###################################################
### code chunk number 5: wordnet.Rnw:91-94 (eval = FALSE)
###################################################
## filter <- getTermFilter("ExactMatchFilter", "company", TRUE)
## terms <- getIndexTerms("NOUN", 1, filter)
## getSynonyms(terms[[1]])


###################################################
### code chunk number 6: wordnet.Rnw:102-103 (eval = FALSE)
###################################################
## synonyms("company", "NOUN")


###################################################
### code chunk number 7: wordnet.Rnw:114-119 (eval = FALSE)
###################################################
## filter <- getTermFilter("ExactMatchFilter", "hot", TRUE)
## terms <- getIndexTerms("ADJECTIVE", 1, filter)
## synsets <- getSynsets(terms[[1]])
## related <- getRelatedSynsets(synsets[[1]], "!")
## sapply(related, getWord)


