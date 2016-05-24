### R code from vignette source 'sos.Rnw'

###################################################
### code chunk number 1: sos.Rnw:10-12
###################################################
options(width = 60, useFancyQuotes = FALSE)
options(repos=c(CRAN="http://cran.cnr.berkeley.edu"))


###################################################
### code chunk number 2: Petal.Length
###################################################
help.search('Petal.Length')


###################################################
### code chunk number 3: PL.RSiteSearch
###################################################
RSiteSearch('Petal.Length')


###################################################
### code chunk number 4: Petal.Length.sos
###################################################
library(sos)
PL <- findFn('Petal.Length')


###################################################
### code chunk number 5: Petal.Length.sos.2
###################################################
PL <- ???Petal.Length


###################################################
### code chunk number 6: summary.PL
###################################################
# the following table has been
# manually edited for clarity
summary(PL)


###################################################
### code chunk number 7: summary.PL-print
###################################################
s <- summary(PL)
blank <- data.frame(Package = "<...>",
                    Count = "", MaxScore = "",
                    TotalScore = "",
                    Date = "")
s$PackageSummary[] <- lapply(s$PackageSummary[], as.character)
row.names(s$PackageSummary) <-
  as.character(s$PackageSummary$Package)
s$PackageSummary <- rbind(s$PackageSummary['yaImpute', ],
                          blank,
                          s$PackageSummary['datasets', ],
                          blank)
print(s, row.names = FALSE)


###################################################
### code chunk number 8: Petal.Length.sos.3
###################################################
PL[PL$Package == 'datasets', 'Function']


###################################################
### code chunk number 9: Petal.Length.sos.3-print
###################################################
print(PL[PL$Package == 'datasets', 'Function'], max.levels = 0)


###################################################
### code chunk number 10: RSiteSearch-spline
###################################################
RSiteSearch('spline')


###################################################
### code chunk number 11: RSiteSearch-spline-numpages
###################################################
getRSiteSearchHits <- function(description) {
  today <- format(Sys.time(), "%Y-%m-%d")
  con <- url(description)
  on.exit(close(con))
  lines <- try(readLines(con))
  if(class(lines) == 'try-error'){
    return(list(hits=0, date=today))
  }
  pattern <- "^.*<!-- HIT -->([0-9]+)<!-- HIT -->.*$"
  hits <- sub(pattern, "\\1", lines[grep(pattern, lines)])
  list(hits = hits, date = today)
}
splineHits <- getRSiteSearchHits("http://search.r-project.org/cgi-bin/namazu.cgi?query=spline&max=20&result=normal&sort=score&idxname=Rhelp08&idxname=functions&idxname=views")


###################################################
### code chunk number 12: RSiteSearch-spline-fun
###################################################
RSiteSearch('spline', 'fun')


###################################################
### code chunk number 13: RSiteSearch-spline-fun-numpages
###################################################
splineFunHits <- getRSiteSearchHits("http://search.r-project.org/cgi-bin/namazu.cgi?query=spline&max=20&result=normal&sort=score&idxname=functions")


###################################################
### code chunk number 14: sos-spline
###################################################
splinePacs <- findFn('spline')


###################################################
### code chunk number 15: sos-spline-maxPages-999
###################################################
splineAll <- findFn('spline', maxPages = 999)


###################################################
### code chunk number 16: sos-spline-subset
###################################################
selSpl <- splineAll[, 'Function'] == 'spline'
splineAll[selSpl, ]


###################################################
### code chunk number 17: sos-spline-grep
###################################################
grepFn('spline', splineAll, ignore.case = TRUE)


###################################################
### code chunk number 18: sos-spline-grep
###################################################
g <- grepFn('spline', splineAll, ignore.case = TRUE)
gFunc6 <- as.character(g[6, "Function"])
gPac6 <- as.character(g[6, "Package"])
gScore6 <- g[6, "Score"]
gCount6 <- g[6, "Count"]


###################################################
### code chunk number 19: writeFindFn2xls-options
###################################################
op <- options(width = 80)


###################################################
### code chunk number 20: writeFindFn2xls
###################################################
writeFindFn2xls(splineAll)


###################################################
### code chunk number 21: writeFindFn2xls-options
###################################################
options(op)


###################################################
### code chunk number 22: install-and-write-options
###################################################
op <- options(width=80)


###################################################
### code chunk number 23: install-and-write
###################################################
splineAll <- findFn('spline', maxPages = 999)
# Do not include in auto test
#installPackages(splineAll)
writeFindFn2xls(splineAll)


###################################################
### code chunk number 24: install-and-write-options-undo
###################################################
options(op)


###################################################
### code chunk number 25: differntial-equations
###################################################
de <- findFn('differential equation')
des <- findFn('differential equations')
de. <- de | des


