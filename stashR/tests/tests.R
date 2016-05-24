## Tests for Regression testing ##########################################
##########################################################################

library(stashR)

wd <- getwd()
dir <- file.path(wd,"testDir")

##########################################################################
## Test objects of class 'remoteDB'
## These tests will fail (within a 'try()') if Internet connectivity
## is not available

myurl <- "http://www.biostat.jhsph.edu/MCAPS/data_v0.3/"

## create a 'remoteDB' object ##
db <- new("remoteDB", url= myurl, dir = dir, name= "MCAPS")
show(db)
show(class(db))
show(db@url)
show(db@dir)


## other prelim steps necessary ##
## dbCreate(db)

## test the methods ##
try( dbList(db) )
x <- try( dbFetch(db, "01073") )
str(x)

try( dbFetch(db, "01004") )
try( dbDelete(db,"01073") )
try( dbInsert(db,key = "01004", value = 1) )

try( dbSync(db) )
dir(file.path(db@dir, "data"))

try( dbSync(db, key = "01073") )
dir(file.path(db@dir, "data"))

try( dbSync(db, key = c("01004","01073")) )
dir(file.path(db@dir, "data"))
try( dbExists(db,c("01073", "01004","55079")) )

## remove db@dir directory ##
unlink(db@dir, recursive = TRUE)

##########################################################################
## Test objects of class 'localDB'

wd <- getwd()
dir <- file.path(wd,"testDir")

## create a 'remoteDB' object ##
dbLocal <- new("localDB", dir= dir, name= "MCAPS")
show(dbLocal)
show(class(dbLocal))
show(dbLocal@dir)

## test the methods  ##
dbInsert(dbLocal,key = "01004", value = 1:10)
dbList(dbLocal)
dbInsert(dbLocal,key = "01005", value = rep(5,10))
dbInsert(dbLocal,key = "01006", value = matrix(1,5,4))
dbList(dbLocal)

reposVersion(dbLocal)
reposVersion(dbLocal) <- 1
dbList(dbLocal)
try( dbFetch(dbLocal, "01005") )
try( dbDelete(dbLocal, "01004") )
try( dbInsert(dbLocal, "01005", 1))

reposVersion(dbLocal) <- -1
dbList(dbLocal)
dbFetch(dbLocal, "01005")

dbFetch(dbLocal, "01004")  
try( dbFetch(dbLocal, "01073") )
dbFetch(dbLocal, "01005")
dbDelete(dbLocal,"01004")
dbList(dbLocal)	
try( dbDelete(dbLocal,"01004") )
dbDelete(dbLocal,"01005")
dbList(dbLocal)
dbExists(dbLocal,key="01004")
dbExists(dbLocal,key="01006")

## Weird object names
dbInsert(dbLocal, "x.1", 1)
dbInsert(dbLocal, "x.1", 2)
dbInsert(dbLocal, "x.2", 3)
dbInsert(dbLocal, "y.1.1.1", 4)

dbList(dbLocal)
dbFetch(dbLocal, "x.2")


dbUnlink(dbLocal)

################################################################################
## Test MD5 digests

db <- new("localDB", dir = "testMD5", name = "testMD5")
dbInsert(db, "obj", rnorm(100))

stopifnot(stashR:::validDataInternal(db, "obj"))

dbUnlink(db)
