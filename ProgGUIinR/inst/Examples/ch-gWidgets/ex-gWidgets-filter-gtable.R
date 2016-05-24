###################################################
### chunk number 1: 
###################################################
#line 18 "ex-gWidgets-filter-gtable.Rnw"
library(gWidgets)


###################################################
### chunk number 2: 
###################################################
#line 25 "ex-gWidgets-filter-gtable.Rnw"
options(repos="http://streaming.stat.iastate.edu/CRAN")


###################################################
### chunk number 3: 
###################################################
#line 28 "ex-gWidgets-filter-gtable.Rnw"
d <- available.packages()       # pick a cran site


###################################################
### chunk number 4: 
###################################################
#line 33 "ex-gWidgets-filter-gtable.Rnw"
w <- gwindow("test of filter")
g <- ggroup(cont=w, horizontal=FALSE)
ed <- gedit("", cont=g)
tbl <- gtable(d, cont=g, filter.FUN="manual", expand=TRUE)


###################################################
### chunk number 5: 
###################################################
#line 49 "ex-gWidgets-filter-gtable.Rnw"
ourMatch <- function(curVal, vals) {
  grepl(curVal, vals)
}


###################################################
### chunk number 6: 
###################################################
#line 61 "ex-gWidgets-filter-gtable.Rnw"
id <- addHandlerKeystroke(ed, handler=function(h, ...) {
  vals <- tbl[, 1, drop=TRUE]
  curVal <- svalue(h$obj)
  vis <- ourMatch(curVal, vals)
  visible(tbl) <- vis
})


