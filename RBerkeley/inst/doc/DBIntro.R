### R code from vignette source 'DBIntro.Rnw'

###################################################
### code chunk number 1: loadRBerkeley
###################################################
library(RBerkeley)
if(file.exists("myDB.db")) {
  unlink("myDB.db")
}


###################################################
### code chunk number 2: dbopen
###################################################
dbh <- db_create()
dbh


###################################################
### code chunk number 3: dbopen
###################################################
ret <- db_open(dbh, txnid=NULL, file="myDB.db", type="BTREE", flags=mkFlags(DB_CREATE, DB_EXCL))
db_strerror(ret)


###################################################
### code chunk number 4: mkFlags
###################################################
mkFlags(DB_CREATE)
mkFlags(DB_EXCL)
mkFlags(DB_CREATE,DB_EXCL)
mkFlags(DB_CREATE,DB_EXCL,DB_EXCL) # bitwise OR duplicates: no change


###################################################
### code chunk number 5: dbput
###################################################
db_put(dbh, key="Ross", data="Ihaka")
db_put(dbh, key="Robert", data="Gentleman")


###################################################
### code chunk number 6: dbputraw
###################################################
db_put(dbh,key=charToRaw("Ross"),data=charToRaw("Ihaka"))
charToRaw("Ihaka")


###################################################
### code chunk number 7: dbcursor
###################################################
dbc <- db_cursor(dbh)


###################################################
### code chunk number 8: dbcursorput
###################################################
dbcursor_put(dbc, key=100L, data=5L, flags=mkFlags(DB_KEYLAST))


###################################################
### code chunk number 9: dbget1
###################################################
db_get(dbh, key=100L)


###################################################
### code chunk number 10: dbget2
###################################################
unserialize(db_get(dbh, key=100L))


###################################################
### code chunk number 11: dbcursorget-1
###################################################
res <- dbcursor_get(dbc, n=1)
res


###################################################
### code chunk number 12: resunserialize
###################################################
lapply(res[[1]], unserialize)


###################################################
### code chunk number 13: dbcursorget-2
###################################################
dbcursor_get(dbc, key="Ross", flags=mkFlags("DB_SET"))


###################################################
### code chunk number 14: dbcursorget-3
###################################################
res <- dbcursor_get(dbc, key="Ross", data="Brawn", flags=mkFlags("DB_SET"))
lapply(res[[1]], unserialize)


###################################################
### code chunk number 15: dbcursorget-4
###################################################
dbcursor_get(dbc, key="Ross", data="Braun", flags=mkFlags("DB_GET_BOTH"))


###################################################
### code chunk number 16: dbcursornew
###################################################
dbcursor_close(dbc)
dbc <- db_cursor(dbh)
res <- dbcursor_get(dbc, flags=mkFlags("DB_NEXT"), n=100)
res


###################################################
### code chunk number 17: dbkeys
###################################################
sapply(res[-1], function(x) unserialize(x$key))


###################################################
### code chunk number 18: dbdel
###################################################
db_strerror(db_exists(dbh, key=charToRaw("Ross")))
db_del(dbh, key=charToRaw("Ross"))
db_strerror(db_exists(dbh, key=charToRaw("Ross")))


###################################################
### code chunk number 19: dbcursordel
###################################################
firstrecord <- dbcursor_get(dbc, flags=mkFlags("DB_FIRST"))[[1]]
db_strerror(dbcursor_del(dbc))
dbcursor_get(dbc, key=firstrecord$key, flags=mkFlags("DB_SET"))


###################################################
### code chunk number 20: close
###################################################
dbcursor_close(dbc)
db_close(dbh)


