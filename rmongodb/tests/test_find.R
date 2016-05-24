library(rmongodb)
library(RUnit)

# 19 tests
# 16.12.2013

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_find", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)

  # inster data
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Dave")
  mongo.bson.buffer.append(buf, "city", "Munich")
  mongo.bson.buffer.append(buf, "age", 27L)
  x <- mongo.bson.from.buffer(buf)
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Fred")
  mongo.bson.buffer.append(buf, "city", "Chicago")
  mongo.bson.buffer.append(buf, "age", 31L)
  y <- mongo.bson.from.buffer(buf)
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Silvia")
  mongo.bson.buffer.append(buf, "city", "London")
  mongo.bson.buffer.append(buf, "age", 24L)
  z <- mongo.bson.from.buffer(buf)
  mongo.insert.batch(mongo, ns, list(x, y, z))
    # insert with JSON
  mongo.insert(mongo, ns, '{"name":"Peter", "city":"Rom", "age":21}')
  
  
  # create bad query
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.start.object(buf, "age")
  mongo.bson.buffer.append(buf, "$bad", 1L)
  mongo.bson.buffer.finish.object(buf)
  query  <- mongo.bson.from.buffer(buf)
  result <- mongo.find.one(mongo, ns, query)
  checkEquals(result, NULL)
  
  
  if (is.null(result)) {
    err <- mongo.get.server.err(mongo)
    print(err)
    err_str <- mongo.get.server.err.string(mongo)
    print(err_str)
  } 
  checkTrue( is.integer(err) )
  checkTrue( is.character(err_str) )
  checkTrue( grep("\\$bad", err_str) == 1 )
  
  
  # good query with find.one
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Dave")
  query  <- mongo.bson.from.buffer(buf)
  result <- mongo.find.one(mongo, ns, query)
  checkEquals( class(result), "mongo.bson")
  checkEquals( mongo.bson.value(result, "name"), "Dave")
  
  # good query with find.one and JSON
  result <- mongo.findOne(mongo, ns, '{"name":"Peter"}')
  checkEquals( class(result), "mongo.bson")
  checkEquals( mongo.bson.value(result, "name"), "Peter")
  
  # good query with find and bson find
  iter <- mongo.bson.find(result, "city")
  if (!is.null(iter)) {
    city <- mongo.bson.iterator.value(iter)
    buf <- mongo.bson.buffer.create()
    mongo.bson.buffer.append(buf, "city", city)
    query <- mongo.bson.from.buffer(buf)
    #print(paste("find city: ", city, sep=""))
    res <- mongo.find.one(mongo, ns, query)
    checkEquals( mongo.bson.value(res, "city"), "Rom")    
  }
  
  
  
  # good query with find and sort
  cursor <- mongo.find(mongo, ns, sort=mongo.bson.from.list(list(city=1L)), limit=100L)
  res <- NULL
  while (mongo.cursor.next(cursor))
    res <- rbind(res, mongo.bson.to.list(mongo.cursor.value(cursor)))
  mongo.cursor.destroy(cursor)  
  checkEquals( class(res), "matrix")
  checkEquals( dim(res), c(4,4))
  checkIdentical( sort(unlist(res[, "city"])) , unlist(res[, "city"]))
  
  # good query with find and sort and JSON
  res <- mongo.find.batch(mongo, ns, '{"age":{"$gt":21}}', sort='{"age":1}', 
                          data.frame=TRUE, mongo.oid2character=TRUE)
  checkEquals( class(res), "data.frame")
  checkTrue( !is.unsorted( res[,"age"] ) )
  
  # good query with find and sort and fields and JSON
  res <- mongo.find.batch(mongo, ns, '{"age":{"$lt":30}}', 
                          sort='{"age":1}', fields='{"city":0}', 
                          data.frame=TRUE, mongo.oid2character=TRUE)
  checkEquals( class(res), "data.frame")
  checkTrue( !is.unsorted( res[,"age"] ) )
  checkEquals( colnames(res), c("_id", "name", "age"))
  
  # good query with find.all
  res <- mongo.find.all(mongo, ns, 
                        data.frame=TRUE, mongo.oid2character=TRUE)  
  checkEquals( dim(res), c(4,4) )
  checkTrue( is.list(res) )
  
  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}