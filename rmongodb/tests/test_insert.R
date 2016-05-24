library(rmongodb)
library(RUnit)

# 7 tests
# 22.11.2013

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_insert", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)
 
  # insert data
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "John")
  mongo.bson.buffer.append(buf, "age", 22L)
  x <- mongo.bson.from.buffer(buf)
  mongo.insert(mongo, ns, x);
  checkEquals(mongo.count(mongo, ns), 1)
  checkEquals( 
    mongo.bson.value(mongo.find.one(mongo, ns, query=list("name"="John")), "name"),
    "John")
  
  
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "John")
  mongo.bson.buffer.append(buf, "age", 27L)
  y <- mongo.bson.from.buffer(buf)
  mongo.insert(mongo, ns, y);  
  checkEquals(mongo.count(mongo, ns), 2)
  checkEquals( 
    mongo.bson.value(mongo.find.one(mongo, ns, query=list("name"="John", "age"=27L)), "age"),
    27)
  
  
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Dwight")
  mongo.bson.buffer.append(buf, "city", "NY")
  z <- mongo.bson.from.buffer(buf)
  mongo.insert(mongo, ns, z)  
  checkEquals(mongo.count(mongo, ns), 3)
  checkEquals( 
    mongo.bson.value(mongo.find.one(mongo, ns, query=list("name"="Dwight")), "name"),
    "Dwight")
  
  
  # insert data with JSON
  mongo.insert(mongo, ns, '{"name":"Peter", "city":{"name":"Rom", "code": 31321}, "age":25}')
  checkEquals(mongo.count(mongo, ns), 4)
  checkEquals( 
    mongo.bson.value(mongo.find.one(mongo, ns, query=list("name"="Peter")), "name"),
    "Peter")
  checkEquals( 
    mongo.bson.value(mongo.find.one(mongo, ns, query=list("name"="Peter")), "city.name"),
    "Rom")
  
  
  # insert batch
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Silvia")
  mongo.bson.buffer.append(buf, "city", "Cincinnati")
  a <- mongo.bson.from.buffer(buf)
  mongo.insert.batch(mongo, ns, list(x, y, z, a))
  checkEquals(mongo.count(mongo, ns), 8)
  
  
  # create index
  mongo.index.create(mongo, ns, "city")
  mongo.index.create(mongo, ns, c("name", "city"))
  # -> more index checks in test_indices

  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}