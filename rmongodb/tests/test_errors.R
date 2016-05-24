library(rmongodb)
library(RUnit)

# 5 tests
# 22.11.2013

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_errors", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)

  mongo.insert(mongo, ns, mongo.bson.from.JSON('{"$or":2}') )
  err <- mongo.get.last.err(mongo, db)
  print(mongo.get.server.err(mongo))
  checkEquals( class(err), "mongo.bson")
  
  iter <- mongo.bson.find(err, "code")
  print(mongo.bson.iterator.value(iter))
  checkTrue( is.integer(iter) )
  
  err_str <- mongo.get.server.err.string(mongo)
  print(err_str)
  checkTrue( is.character(err_str) )
  
  iter <- mongo.bson.find(err, "err")
  err_str <- mongo.bson.iterator.value(iter)
  print(err_str)
  checkTrue( is.character(err_str) )
  
  out <- mongo.reset.err(mongo, db)
  checkTrue( is.null( out ) )
  
  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}