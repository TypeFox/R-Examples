## ----connect-------------------------------------------------------------
library(rmongodb)
mongo <- mongo.create()
mongo.is.connected(mongo)

## ----insertZips----------------------------------------------------------
# load example data set from rmongodb
data(zips)
head(zips)
zips[1,]$loc

# rename _id field. The original zips data set holds duplicate _id values which will fale during the import
colnames(zips)[5] <- "orig_id"

# create BSON batch object
ziplist <- list()
ziplist <- apply( zips, 1, function(x) c( ziplist, x ) )
res <- lapply( ziplist, function(x) mongo.bson.from.list(x) )

if(mongo.is.connected(mongo) == TRUE){
  mongo.insert.batch(mongo, "rmongodb.zips", res )
}

## ----insertCheck---------------------------------------------------------
dim(zips)
if(mongo.is.connected(mongo) == TRUE){
  nr <- mongo.count(mongo, "rmongodb.zips")
  print( nr )
  res <- mongo.find.all(mongo, "rmongodb.zips", limit=20)
  head( res )
}

## ----aggregationFramework------------------------------------------------
pipe_1 <- mongo.bson.from.JSON('{"$group":{"_id":"$state", "totalPop":{"$sum":"$pop"}}}')
cmd_list <- list(pipe_1)
cmd_list
if(mongo.is.connected(mongo) == TRUE){
  res <- mongo.aggregation(mongo, "rmongodb.zips", cmd_list)
  head( mongo.bson.value(res, "result") )
}

## ----aggregationFramework2-----------------------------------------------
pipe_1 <- mongo.bson.from.JSON('{"$group":{"_id":"$state", "totalPop":{"$sum":"$pop"}}}')
pipe_2 <- mongo.bson.from.JSON('{"$match":{"totalPop":{"$gte":15000000}}}')
cmd_list <- list(pipe_1, pipe_2)
if(mongo.is.connected(mongo) == TRUE){
  res <- mongo.aggregation(mongo, "rmongodb.zips", cmd_list)
  res
}

## ----gridfs--------------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE){
  mgrids <- mongo.gridfs.create(mongo, "rmongodb", prefix = "fs")
  mongo.gridfs.store.file(mgrids, "faust.txt", "Faust")
  gf <- mongo.gridfs.find(mgrids, "Faust")
  mongo.gridfile.get.length(gf)
  mongo.gridfile.get.chunk.count(gf)
}

## ----drop----------------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE){
  mongo.drop(mongo, "rmongodb.zips")
  mongo.drop.database(mongo, "rmongodb")
  
  # close connection
  mongo.destroy(mongo)
}

