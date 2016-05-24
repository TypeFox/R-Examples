library(rmongodb)
library(RUnit)

# 2 tests
# 22.11.2013


# test create bson object with values
r <- as.integer(c(1,2,3,4,5,6,7,8))
dim(r) <- c(2,2,2)
t <- Sys.time()
t <- c(t, t+10, t+60, t+120)
dim(t) <- c(2,2)
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append.int(buf, "test", r)
mongo.bson.buffer.append(buf, "times", t)
b <- mongo.bson.from.buffer(buf)

# iterate over bson object
iter <- mongo.bson.iterator.create(b)
while (mongo.bson.iterator.next(iter)) {
  print(mongo.bson.iterator.key(iter))
  print(mongo.bson.iterator.value(iter))
}

iter <- mongo.bson.find(b, "test")
out <- mongo.bson.iterator.value(iter)
checkEquals( out, r )

sub2 <- mongo.bson.value(b, "times")
checkEqualsNumeric( sub2, t, tolerance=1.0e-6)
