library(filehash)

name <- sprintf("testdb-v%s", packageDescription("filehash", fields = "Version"))
dbCreate(name, "DB1")
db <- dbInit(name, "DB1")

set.seed(1)
dbInsert(db, "a", rnorm(10))
dbInsert(db, "b", runif(7))
dbInsert(db, "list", list(1, 2, 3, 4, 5, 6, "a"))
dbInsert(db, "c", 1L)
dbInsert(db, "entry", "string")
dbDelete(db, "b")
         
