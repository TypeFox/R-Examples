library(rmongodb)
library(RUnit)

# 3 tests
# 22.11.2013

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_indices", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)
   
#   for( i in rep(1:10,3)){
#     mongo.insert(mongo, ns, mongo.bson.from.JSON(paste('{"b":',i,'}')))
#   }
#   count1 <- mongo.count(mongo, ns)
#   out <- mongo.index.create(mongo, ns, "b", options=mongo.index.unique)
#   checkTrue( is.null(out) )
#   count2 <- mongo.count(mongo, ns)
#   checkTrue( count2 < count1)
  
  print("test add dup key")
  for( i in 1:10){
    mongo.insert(mongo, ns, mongo.bson.from.JSON(paste('{"a":',i,', "b":"', letters[i],'"}', sep="")))
  }
  out <- mongo.index.create(mongo, ns, "a", mongo.index.unique)
  checkTrue( is.null(out) )
  insert <- mongo.insert(mongo, ns, mongo.bson.from.JSON('{"a":5}'))
  checkTrue( !insert )
  
  # create indey by key
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "b", "")
  mongo.bson.buffer.append(buf, "a", 1L)
  b <- mongo.bson.from.buffer(buf)
  print("index create")
  out  <- mongo.index.create(mongo, ns, b)
  checkTrue( is.null(out) )
  
  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}