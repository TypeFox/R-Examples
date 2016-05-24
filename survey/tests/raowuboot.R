## regression test for bug reported by Richard Valliant
library(survey)
s<-subbootweights(c(1,1),1:2, 50)
stopifnot(all(s$repweights$weights %in% c(0,2)))
