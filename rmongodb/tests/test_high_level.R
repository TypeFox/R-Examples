library(rmongodb)
library(RUnit)

# 03.03.2014

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_high_level", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)
  # and insert some data
  mongo.insert(mongo, ns, '{"name":"Peter", "city":"Rom"}')
  mongo.insert(mongo, ns, '{"name":"Markus", "city":"Munich", "age":51}')
  mongo.insert(mongo, ns, '{"name":"Tom", "city":"London", "age":1}')
  mongo.insert(mongo, ns, '{"name":"Jon", "age":23}')
  
  
  #mongo.distinct
  res <- mongo.distinct(mongo, ns, "name")
  checkEquals(as.vector(res), c("Peter", "Markus", "Tom", "Jon"))
  
  res <- mongo.distinct(mongo, ns, "name2")
  checkEquals(length(res), 0)
    
  res <- mongo.distinct(mongo, ns, "age")
  checkEquals(as.vector(res), c(51,1,23))
  
  checkException( mongo.distinct(mongo, "ns", "age"), "Wrong namespace (ns)." )
  
  
  # mongo.parse.ns
  res <- rmongodb:::mongo.parse.ns("db.coll")
  checkEquals(res$db, "db")
  checkEquals(res$collection, "coll")
  checkTrue(is.list(res))
  
  res <- rmongodb:::mongo.parse.ns("db_coll")
  checkTrue(is.null(res))
  
  
  #mongo.aggregation
  # load more data for good tests
  mongo.drop(mongo, ns)
  data(zips)
  colnames(zips)[5] <- "orig_id"
  ziplist <- list()
  ziplist <- apply( zips, 1, function(x) c( ziplist, x ) )
  res <- lapply( ziplist, function(x) mongo.bson.from.list(x) )
  mongo.insert.batch(mongo, ns, res )
  
  pipe_1 <- mongo.bson.from.JSON('{"$group":{"_id":"$state", "totalPop":{"$sum":"$pop"}}}')
  cmd_list <- list(pipe_1)
  res <- mongo.aggregation(mongo, ns, cmd_list)
  res <- mongo.bson.value(res, "result")
  checkEquals(length(res), 51)
  
  checkException( mongo.aggregation(mongo, "ns", cmd_list), "Wrong namespace (ns)." )
  
  pipe_1 <- mongo.bson.from.JSON('{"$group":{"_ids":"$state", "totalPop":{"$sum":"$pop"}}}')
  cmd_list <- list(pipe_1)
  checkException( mongo.aggregation(mongo, ns, cmd_list), "mongoDB error: 10. Please check ?mongo.get.err for more details." )
  
  pipe_1 <- mongo.bson.from.JSON('{"$group":{"_id":"$state", "totalPop":{"$sum":"$pop"}}}')
  pipe_2 <- mongo.bson.from.JSON('{"$match":{"totalPop":{"$gte":10000000}}}')
  cmd_list <- list(pipe_1, pipe_2)
  res <- mongo.aggregation(mongo, ns, cmd_list)
  res <- mongo.bson.value(res, "result")
  checkEquals(length(res), 7)
  
  
  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}

