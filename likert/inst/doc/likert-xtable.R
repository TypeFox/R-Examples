### R code from vignette source 'likert-xtable.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(digits=3)
options(width=80)
options(continue="  ")


###################################################
### code chunk number 2: likert-xtable.Rnw:22-42
###################################################
library(likert)
library(ggplot2)
data(pisaitems)

##### Item 24: Reading Attitudes
items24 <- pisaitems[,substr(names(pisaitems), 1,5) == 'ST24Q']
names(items24) <- c(
	  "I read only if I have to.",
	  "Reading is one of my favorite hobbies.",
	  "I like talking about books with other people.",
	  "I find it hard to finish books.",
	  "I feel happy if I receive a book as a present.",
	  "For me, reading is a waste of time.",
	  "I enjoy going to a bookstore or a library.",
	  "I read only to get information that I need.",
	  "I cannot sit still and read for more than a few minutes.",
	  "I like to express my opinions about books I have read.",
	  "I like to exchange books with my friends.")
l24 = likert(items24)
l24g <- likert(items24, grouping=pisaitems$CNT)


###################################################
### code chunk number 3: likert-xtable.Rnw:47-48
###################################################
xtable(l24)


###################################################
### code chunk number 4: likert-xtable.Rnw:51-52 (eval = FALSE)
###################################################
## xtable(l24g)


