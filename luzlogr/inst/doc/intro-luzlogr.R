## ------------------------------------------------------------------------
library(luzlogr)
openlog("test.log")
printlog("message")
closelog()

## ---- echo=FALSE---------------------------------------------------------
lg <- readLines("test.log")
invisible(file.remove("test.log"))
for(i in seq_along(lg)) cat(lg[i], "\n")

## ------------------------------------------------------------------------
openlog("test.log")
printlog("message", 1, 2)
printlog(head(cars))
closelog(sessionInfo = FALSE)

## ---- echo=FALSE---------------------------------------------------------
lg <- readLines("test.log")
invisible(file.remove("test.log"))
for(i in seq_along(lg)) cat(lg[i], "\n")

## ------------------------------------------------------------------------
openlog("test.log", loglevel = 0)
printlog("This message will appear", level = 0)
printlog("So will this (level 0 by default)")
printlog("This will not", level = -1)
closelog(sessionInfo = FALSE)

## ---- echo=FALSE---------------------------------------------------------
lg <- readLines("test.log")
invisible(file.remove("test.log"))
for(i in seq_along(lg)) cat(lg[i], "\n")

## ------------------------------------------------------------------------
openlog("test.log")
printlog("A normal message")
printlog("A flagged message!", flag = TRUE)
flaglog("Another")
closelog(sessionInfo = FALSE)

## ---- echo=FALSE---------------------------------------------------------
lg <- readLines("test.log")
invisible(file.remove("test.log"))
for(i in seq_along(lg)) cat(lg[i], "\n")

## ------------------------------------------------------------------------
con <- gzfile("test.log.gz")
openlog(con)
printlog("Sending to a compressed logfile")
closelog(sessionInfo = FALSE)

## ---- echo=FALSE---------------------------------------------------------
invisible(file.remove("test.log.gz"))

