
### ceeboo 2012

ns <- "arulesSequences"
library(ns, character.only = TRUE)

k <- objects(envir = asNamespace(ns), pattern = "^R_")
names(k) <- k
k <- sapply(lapply(k, get, envir = asNamespace(ns)), inherits, 
		   "NativeSymbolInfo")
stopifnot(all(k))

###
