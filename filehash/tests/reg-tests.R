suppressMessages(library(filehash))

######################################################################
## Test 'filehashRDS' class

dbCreate("mydbRDS", "RDS")
db <- dbInit("mydbRDS", "RDS")
show(db)

## Put some data into it
set.seed(1000)
dbInsert(db, "a", 1:10)
dbInsert(db, "b", rnorm(100))
dbInsert(db, "c", 100:1)
dbInsert(db, "d", runif(1000))
dbInsert(db, "other", "hello")

dbList(db)

dbExists(db, "e")
dbExists(db, "a")

env <- db2env(db)
ls(env)

env$a
env$b
env$c
str(env$d)
env$other

env$b <- rnorm(100)
mean(env$b)

env$a[1:5] <- 5:1
print(env$a)

dbDelete(db, "c")

tryCatch(print(env$c), error = function(e) cat(as.character(e)))
tryCatch(dbFetch(db, "c"), error = function(e) cat(as.character(e)))

## Check trailing '/' problem
dbCreate("testRDSdb", "RDS")
db <- dbInit("testRDSdb/", "RDS")
print(db)

######################################################################
## test filehashDB1 class

dbCreate("mydb", "DB1")
db <- dbInit("mydb", "DB1")

## Put some data into it
set.seed(1000)
dbInsert(db, "a", 1:10)
dbInsert(db, "b", rnorm(100))
dbInsert(db, "c", 100:1)
dbInsert(db, "d", runif(1000))
dbInsert(db, "other", "hello")

dbList(db)

env <- db2env(db)
ls(env)

env$a
env$b
env$c
str(env$d)
env$other

env$b <- rnorm(100)
mean(env$b)

env$a[1:5] <- 5:1
print(env$a)

dbDelete(db, "c")

tryCatch(print(env$c), error = function(e) cat(as.character(e)))
tryCatch(dbFetch(db, "c"), error = function(e) cat(as.character(e)))

numbers <- rnorm(100)
dbInsert(db, "numbers", numbers)
b <- dbFetch(db, "numbers")
stopifnot(all.equal(numbers, b))
stopifnot(identical(numbers, b))

################################################################################
## Other tests

rm(list = ls())


dbCreate("testLoadingDB", "DB1")
db <- dbInit("testLoadingDB", "DB1")

set.seed(234)

db$a <- rnorm(100)
db$b <- runif(1000)

dbLoad(db)  ## 'a', 'b'
summary(a)
summary(b)

rm(list = ls())
db <- dbInit("testLoadingDB", "DB1")

dbLazyLoad(db)

summary(a)
summary(b)



################################################################################
## Check dbReorganize

dbCreate("test_reorg", "DB1")
db <- dbInit("test_reorg", "DB1")

set.seed(1000)
dbInsert(db, "a", 1)
dbInsert(db, "a", 1)
dbInsert(db, "a", 1)
dbInsert(db, "a", 1)
dbInsert(db, "b", rnorm(1000))
dbInsert(db, "b", rnorm(1000))
dbInsert(db, "b", rnorm(1000))
dbInsert(db, "b", rnorm(1000))
dbInsert(db, "c", runif(1000))
dbInsert(db, "c", runif(1000))
dbInsert(db, "c", runif(1000))
dbInsert(db, "c", runif(1000))

summary(db$b)
summary(db$c)

print(file.info(db@datafile)$size)

dbReorganize(db)

db <- dbInit("test_reorg", "DB1")

print(file.info(db@datafile)$size)

summary(db$b)
summary(db$c)


################################################################################
## Taken from the vignette

file.remove("mydb")

dbCreate("mydb")
db <- dbInit("mydb")

set.seed(100)

dbInsert(db, "a", rnorm(100))
value <- dbFetch(db, "a")
mean(value)

dbInsert(db, "b", 123)
dbDelete(db, "a")
dbList(db)
dbExists(db, "a")

file.remove("mydb")

################################################################################
## Check queue

db <- createQ("testq")
push(db, 1)
push(db, 2)
top(db)

pop(db)
top(db)
