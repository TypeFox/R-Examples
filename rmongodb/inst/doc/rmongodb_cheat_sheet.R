### R code from vignette source 'rmongodb_cheat_sheet.Rnw'

###################################################
### code chunk number 1: rmongodb_cheat_sheet.Rnw:81-82 (eval = FALSE)
###################################################
## install.packages("rmongodb")


###################################################
### code chunk number 2: rmongodb_cheat_sheet.Rnw:85-87 (eval = FALSE)
###################################################
## library(devtools)
## install_github("rmongodb", "mongosoup")


###################################################
### code chunk number 3: rmongodb_cheat_sheet.Rnw:90-91 (eval = FALSE)
###################################################
## library(rmongodb)


###################################################
### code chunk number 4: rmongodb_cheat_sheet.Rnw:94-95 (eval = FALSE)
###################################################
## ??rmongodb


###################################################
### code chunk number 5: rmongodb_cheat_sheet.Rnw:103-104 (eval = FALSE)
###################################################
## mongo <- mongo.create()


###################################################
### code chunk number 6: rmongodb_cheat_sheet.Rnw:107-110 (eval = FALSE)
###################################################
## mongo <- mongo.create(host="127.1.1.1:27017", 
##             username="USER", password="XXX", 
##             db="database")


###################################################
### code chunk number 7: rmongodb_cheat_sheet.Rnw:113-114 (eval = FALSE)
###################################################
## mongo.is.connected(mongo)


###################################################
### code chunk number 8: rmongodb_cheat_sheet.Rnw:117-118 (eval = FALSE)
###################################################
## mongo.destroy(mongo)


###################################################
### code chunk number 9: rmongodb_cheat_sheet.Rnw:125-135 (eval = FALSE)
###################################################
## data(zips)
## # rename _id field. The original zips data set holds duplicate _id values which will fale during the import
## colnames(zips)[5] <- "orig_id"
## # create BSON batch object
## ziplist <- list()
## ziplist <- apply( zips, 1, function(x) c( ziplist, x ) )
## res <- lapply( ziplist, function(x) mongo.bson.from.list(x) )
## if(mongo.is.connected(mongo) == TRUE){
##   mongo.insert.batch(mongo, "rmongodb.zips", res )
## }


###################################################
### code chunk number 10: rmongodb_cheat_sheet.Rnw:143-145 (eval = FALSE)
###################################################
## mongo.get.databases(mongo)
## mongo.get.database.collections(mongo, "rmongodb")


###################################################
### code chunk number 11: rmongodb_cheat_sheet.Rnw:148-151 (eval = FALSE)
###################################################
## mongo.get.err(mongo)
## mongo.get.server.err(mongo)
## mongo.get.server.err.string(mongo)


###################################################
### code chunk number 12: rmongodb_cheat_sheet.Rnw:154-156 (eval = FALSE)
###################################################
## mongo.get.primary(mongo)
## mongo.get.hosts(mongo)


###################################################
### code chunk number 13: rmongodb_cheat_sheet.Rnw:164-165 (eval = FALSE)
###################################################
## mongo.count(mongo, "rmongodb.zips")


###################################################
### code chunk number 14: rmongodb_cheat_sheet.Rnw:170-171 (eval = FALSE)
###################################################
## mongo.get.values(mongo, "rmongodb.zips", "city")


###################################################
### code chunk number 15: rmongodb_cheat_sheet.Rnw:179-181 (eval = FALSE)
###################################################
## mongo.count(mongo, "rmongodb.zips", 
##             query='{"state":"AL"}')


###################################################
### code chunk number 16: rmongodb_cheat_sheet.Rnw:185-187 (eval = FALSE)
###################################################
## bson <- mongo.findOne(mongo, "rmongodb.zips", 
##                       query='{"state":"AL"}')


###################################################
### code chunk number 17: rmongodb_cheat_sheet.Rnw:190-192 (eval = FALSE)
###################################################
## cursor <- mongo.find(mongo, "rmongodb.zips", 
##                      query='{"state":"AL"}')


###################################################
### code chunk number 18: rmongodb_cheat_sheet.Rnw:195-196 (eval = FALSE)
###################################################
## mongo.cursor.to.list(cursor)


###################################################
### code chunk number 19: rmongodb_cheat_sheet.Rnw:199-201 (eval = FALSE)
###################################################
## mongo.find.all(mongo, "rmongodb.zips", 
##               query='{"state":"AL"}')


###################################################
### code chunk number 20: rmongodb_cheat_sheet.Rnw:207-210 (eval = FALSE)
###################################################
## mongo.find.all(mongo, "rmongodb.zips", 
##               query='{"state":"AL"}',
##               skip=5, limit=10)


###################################################
### code chunk number 21: rmongodb_cheat_sheet.Rnw:213-217 (eval = FALSE)
###################################################
## mongo.find.all(mongo, "rmongodb.zips", 
##               query='{"state":"AL"}',
##               fields='{"city":1, "pop":1, "_id":0}', 
##               sort='{"pop":1}')


###################################################
### code chunk number 22: rmongodb_cheat_sheet.Rnw:223-225 (eval = FALSE)
###################################################
## mongo.find.all(mongo, "rmongodb.zips", 
##               query='{"state":"AL", "city":"ACMAR"}')


###################################################
### code chunk number 23: rmongodb_cheat_sheet.Rnw:229-231 (eval = FALSE)
###################################################
## mongo.find.all(mongo, "rmongodb.zips", 
##               query='{"pop":{"$gte":80000}}')


###################################################
### code chunk number 24: rmongodb_cheat_sheet.Rnw:235-237 (eval = FALSE)
###################################################
## mongo.count(mongo, "rmongodb.zips", 
##               query='{"loc":{"$exists":1}}')


###################################################
### code chunk number 25: rmongodb_cheat_sheet.Rnw:249-250 (eval = FALSE)
###################################################
## mongo.bson.to.list(bson)


###################################################
### code chunk number 26: rmongodb_cheat_sheet.Rnw:253-254 (eval = FALSE)
###################################################
## mongo.bson.value(bson, "state")


###################################################
### code chunk number 27: rmongodb_cheat_sheet.Rnw:260-261 (eval = FALSE)
###################################################
## mongo.bson.from.JSON('{"state":"AL"}')


###################################################
### code chunk number 28: rmongodb_cheat_sheet.Rnw:264-267 (eval = FALSE)
###################################################
## buf <- mongo.bson.buffer.create()
## mongo.bson.buffer.append(buf, "state", "AL")
## b <- mongo.bson.from.buffer(buf)


###################################################
### code chunk number 29: rmongodb_cheat_sheet.Rnw:270-271 (eval = FALSE)
###################################################
## ?mongo.bson


###################################################
### code chunk number 30: rmongodb_cheat_sheet.Rnw:279-281 (eval = FALSE)
###################################################
## mongo.insert(mongo, "rmongodb.insert", 
##              '{"user":"markus", "city":"munich"}')


###################################################
### code chunk number 31: rmongodb_cheat_sheet.Rnw:284-290 (eval = FALSE)
###################################################
## bson1 <- mongo.bson.from.JSON(
##   '{"user":"markus", "city":"munich"}')
## bson2 <- mongo.bson.from.JSON(
##   '{"user":"peter", "city":"New York"}')
## mongo.insert.batch(mongo, "rmongodb.insert",
##                    list(bson1, bson2))


###################################################
### code chunk number 32: rmongodb_cheat_sheet.Rnw:293-294 (eval = FALSE)
###################################################
## mongo.index.create(mongo, "rmongodb.insert", '{"user":1}')


###################################################
### code chunk number 33: rmongodb_cheat_sheet.Rnw:297-300 (eval = FALSE)
###################################################
## mongo.update(mongo, "rmongodb.insert", 
##              '{"user":"markus"}', 
##              '{"user":"markus", "city":"berlin"}')


###################################################
### code chunk number 34: rmongodb_cheat_sheet.Rnw:307-315 (eval = FALSE)
###################################################
## pipe_1 <- mongo.bson.from.JSON('{"$group":
##                                {"_id":"$state", "totalPop":
##                                {"$sum":"$pop"}}}')
## pipe_2 <- mongo.bson.from.JSON('{"$match":
##                                {"totalPop":
##                                {"$gte":15000000}}}')
## cmd_list <- list(pipe_1, pipe_2)
## bson <- mongo.aggregation(mongo, "rmongodb.zips", cmd_list)


###################################################
### code chunk number 35: rmongodb_cheat_sheet.Rnw:323-325 (eval = FALSE)
###################################################
## mongo.drop.database(mongo, "rmongodb")
## mongo.destroy(mongo)


