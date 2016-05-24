
### ceeboo 2014

library("arulesSequences")

##
z <- list(a = c("A","B"), b = c("D"), c = c("A","H"))
z <- as(z, "tidLists")
summary(z)
transactionInfo(z)

###
