## ----installRmongodb, eval=FALSE-----------------------------------------
#  install.packages("rmongodb")

## ----installDEV, eval=FALSE----------------------------------------------
#  library(devtools)
#  install_github(repo = "mongosoup/rmongodb")

## ----loadRmongodb--------------------------------------------------------
library(rmongodb)

## ----connect2Mongo-------------------------------------------------------
help("mongo.create")
mongo <- mongo.create()
mongo
mongo.is.connected(mongo)

## ----importZIPdata, echo=FALSE, warning=FALSE, results='hide'------------
if(mongo.is.connected(mongo) == TRUE) {
  # load some data
  library(jsonlite)
  data(zips)
  # rename _id field. The original zips data set holds duplicate _id values which will fale during the import
  colnames(zips)[5] <- "orig_id"
  
  ziplist <- list()
  ziplist <- apply( zips, 1, function(x) c( ziplist, x ) )
  res <- lapply( ziplist, function(x) mongo.bson.from.list(x) )
  
  mongo.insert.batch(mongo, "rmongodb.zips", res )
}

## ----getDBs--------------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  mongo.get.databases(mongo)
}

## ----getColls------------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  db <- "rmongodb"
  mongo.get.database.collections(mongo, db)
}
coll <- "rmongodb.zips"

## ----count, echo=TRUE----------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  help("mongo.count")
  mongo.count(mongo, coll)
}

## ----findOneFirst, echo=TRUE---------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  mongo.find.one(mongo, coll)
}

## ----Distinct, echo=TRUE-------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  res <- mongo.distinct(mongo, coll, "city")
  head(res, 2)
}

## ----findOne-------------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  cityone <- mongo.find.one(mongo, coll, '{"city":"COLORADO CITY"}')
  print( cityone )
  mongo.bson.to.list(cityone)
}

## ----createBSONfromList--------------------------------------------------
query <- mongo.bson.from.list(list('city' = 'COLORADO CITY'))
query
query <- mongo.bson.from.list(list('city' = 'COLORADO CITY', 'loc' = list(-112.952427, 36.976266)))
query

## ----createBSON----------------------------------------------------------
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append(buf, "city", "COLORADO CITY")
query <- mongo.bson.from.buffer(buf)
query

## ----createBSONoneLine---------------------------------------------------
mongo.bson.from.JSON('{"city":"COLORADO CITY", "loc":[-112.952427, 36.976266]}')

## ----createBSONwithDate--------------------------------------------------
date_string <- "2014-10-11 12:01:06"
# Pay attention to timezone argument
query <- mongo.bson.from.list(list(date = as.POSIXct(date_string, tz='MSK')))
# Note, that internall MongoDB strores dates in unixtime format:
query

## ----findMore, warning=FALSE---------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  pop <- mongo.distinct(mongo, coll, "pop")
  hist(pop)
  boxplot(pop)

  nr <- mongo.count(mongo, coll, list('pop' = list('$lte' = 2)))
  print( nr )
  pops <- mongo.find.all(mongo, coll, list('pop' = list('$lte' = 2)))
  head(pops, 2)
}

## ----complexQuery--------------------------------------------------------
# or do it R way, as recommended above:
if(mongo.is.connected(mongo) == TRUE) {
  pops1 <- mongo.find.all(mongo, coll, query = list('pop' = list('$lte' = 2), 'pop' = list('$gte' = 1)))
  head(pops1, 2)
}


## ----complexQueryJSON----------------------------------------------------
library(jsonlite)
json <- '{"pop":{"$lte":2}, "pop":{"$gte":1}}'
cat(prettify(json))
validate(json)
if(mongo.is.connected(mongo) == TRUE) {
  pops1 <- mongo.find.all(mongo, coll, query = list('pop' = list('$lte' = 2), 'pop' = list('$gte' = 1)))
  pops2 <- mongo.find.all(mongo, coll, json)
  identical(pops1, pops2)
}

## ----insert--------------------------------------------------------------
# insert data
a <- mongo.bson.from.JSON( '{"ident":"a", "name":"Markus", "age":33}' )
b <- mongo.bson.from.JSON( '{"ident":"b", "name":"MongoSoup", "age":1}' )
c <- mongo.bson.from.JSON( '{"ident":"c", "name":"UseR", "age":18}' )

if(mongo.is.connected(mongo) == TRUE) {
  icoll <- paste(db, "test", sep=".")
  mongo.insert.batch(mongo, icoll, list(a,b,c) )

  dbs <- mongo.get.database.collections(mongo, db)
  print(dbs)
  mongo.find.all(mongo, icoll)
}

## ----update--------------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  mongo.update(mongo, icoll, list('ident' = 'b'), list('$inc' = list('age' = 3)))

  res <- mongo.find.all(mongo, icoll)
  print(res)
  
  # Creating an index for the field 'ident'
  mongo.index.create(mongo, icoll, list('ident' = 1))
  # check mongoshell!
}

## ----dropColls-----------------------------------------------------------
if(mongo.is.connected(mongo) == TRUE) {
  mongo.drop(mongo, icoll)
  #mongo.drop.database(mongo, db)
  res <- mongo.get.database.collections(mongo, db)
  print(res)
  
  # close connection
  mongo.destroy(mongo)
}

