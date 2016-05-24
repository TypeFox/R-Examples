suppressMessages(library(filehashSQLite))

dbCreate("test1", "SQLite")
db <- dbInit("test1", "SQLite")

set.seed(234)
val <- rnorm(1000)

dbInsert(db, "a", val)
x <- dbFetch(db, "a")

str(x)

stopifnot(identical(x, val))

dbDisconnect(db)
