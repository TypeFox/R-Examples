library(rmongodb)
library(RUnit)

# 1 tests
# 22.11.2013

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_update", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)
  
  # inster data
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Dave")
  mongo.bson.buffer.append(buf, "age", 27L)
  x <- mongo.bson.from.buffer(buf)
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Fred")
  mongo.bson.buffer.append(buf, "age", 31L)
  y <- mongo.bson.from.buffer(buf)
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Silvia")
  mongo.bson.buffer.append(buf, "age", 24L)
  z <- mongo.bson.from.buffer(buf)
  mongo.insert.batch(mongo, ns, list(x, y, z))
  
  
  # update one document
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Silvia")
  query <- mongo.bson.from.buffer(buf)
  
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.start.object(buf, "$inc")
  mongo.bson.buffer.append(buf, "age", 1L)
  mongo.bson.buffer.finish.object(buf)
  op  <- mongo.bson.from.buffer(buf)
  
  mongo.update(mongo, ns, query, op)
   
  checkEquals(
      mongo.bson.value( mongo.find.one(mongo, ns, query), "age"),
      mongo.bson.value( z, "age") +1 )
  
  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}