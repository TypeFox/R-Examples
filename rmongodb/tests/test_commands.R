library(rmongodb)
library(RUnit)

# 18 tests
# 22.11.2013

# set up mongoDB connection and db / collection parameters
mongo <- mongo.create()
db <- "rmongodb"
ns <- paste(db, "test_commands", sep=".")

if( mongo.is.connected(mongo) ){
  
  # clean up old existing collection
  mongo.drop(mongo, ns)
  
  print(mongo.get.primary(mongo))
  test <- mongo.is.master(mongo)
  print(sprintf("IsMaster (%s)", if (test) "Yes" else "No"))
  checkTrue(test)
  
  t <- 2000
  mongo.set.timeout(mongo, t)
  tout <- mongo.get.timeout(mongo)
  checkEquals(t, tout)
  
  out <- mongo.simple.command(mongo, "admin", "buildInfo", 1)
  print(out)
  checkEquals( class(out), "mongo.bson")
  checkEquals( mongo.bson.value(out, "ok"), 1 )
  
  out <- mongo.get.databases(mongo)
  print(out)
  checkTrue( any(out=="rmongodb") )
  
  out <- mongo.simple.command(mongo, "admin", "top", 1L) 
  checkEquals( mongo.bson.value(out, "ok"), 1 )
  checkEquals( class(out), "mongo.bson")
  
  checkTrue(!mongo.drop(mongo, ns))
  mongo.insert(mongo, ns, mongo.bson.from.JSON('{"a":5}'))
  print(paste("drop collection", ns))
  checkTrue(mongo.drop(mongo, ns))

  print("drop database")
  checkTrue(mongo.drop.database(mongo, db))
  
  print("add user")
  checkTrue(mongo.add.user(mongo, "Gerald", "PaSsWoRd"))
  
  # create error
  mongo.simple.command(mongo, db, "badcommand", 0L)
  err <- mongo.get.err(mongo)
  print(err)
  checkTrue( is.integer(err) )
  mongo.reset.err(mongo, db)
  
  print("rename collection")
  for( i in 1:100){
    mongo.insert(mongo, ns, mongo.bson.from.JSON(paste('{"a":',i,'}')))
  }
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "renameCollection", ns)
  mongo.bson.buffer.append(buf, "to", "rmongodb.test_new")
  command <- mongo.bson.from.buffer(buf)
  out <- mongo.command(mongo, "admin", command)
  print(out)
  checkEquals( class(out), "mongo.bson")
  checkEquals( mongo.bson.value(out, "ok"), 1 )
  
  ns <- "rmongodb.test_new"
  
  print("count check")
  count1 <- mongo.count(mongo, ns)
  print(count1)
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "count", "test_new")
  mongo.bson.buffer.append(buf, "query", mongo.bson.empty())
  command <- mongo.bson.from.buffer(buf)
  count2 <- mongo.command(mongo, "rmongodb", command)
  print(count2)
  checkEquals( class(count2), "mongo.bson")
  checkEquals( mongo.bson.value(count2, "ok"), 1 )
  checkEquals( count1, mongo.bson.value(count2, "n"))
  
  
  buf <- mongo.bson.buffer.create()
  mongo.bson.buffer.append(buf, "name", "Ford")
  mongo.bson.buffer.append(buf, "engine", "Vv8")
  z <- mongo.bson.from.buffer(buf)
  mongo.insert(mongo, "rmongodb.test_new2", z)
  out <- mongo.get.database.collections(mongo, "rmongodb")
  checkEquals( length(out), 2)
  
  # cleanup db and close connection
  mongo.drop.database(mongo, db)
  mongo.destroy(mongo)
}