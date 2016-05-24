library(rmongodb)
library(RUnit)
library(jsonlite)

# 24.02.2014

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
checkEquals(mongo.bson.value(b, "test"), r)
checkEqualsNumeric(mongo.bson.value(b, "times"), t, tolerance=1.0e-6)


# test for issue #18
tmp <- mongo.bson.from.list(list(a=TRUE))
checkTrue( mongo.bson.value(tmp, "a") )
print(mongo.bson.to.list(tmp))

tmp <- mongo.bson.from.list(list(a=FALSE))
checkTrue( !mongo.bson.value(tmp, "a") )
print( mongo.bson.to.list(tmp) )

buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append(buf, "test", TRUE)
b <- mongo.bson.from.buffer(buf)
checkTrue(mongo.bson.value(b, "test"))

buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append(buf, "test", FALSE)
b <- mongo.bson.from.buffer(buf)
checkTrue(!mongo.bson.value(b, "test"))


# more test for bson.to.list and to.Robject
b <- mongo.bson.from.JSON('{"First":"Joe Smith", "Last":21.5}')
l <- mongo.bson.to.Robject(b)
checkTrue(is.list(l))

b <- mongo.bson.from.JSON('{"First":"Joe", "Last":"Smith"}')
l <- mongo.bson.to.Robject(b)
checkTrue(is.vector(l))


# test create bson oject with raw data
r <- as.raw(r)
dim(r) <- c(2,4)
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append(buf, "test", r)
b <- mongo.bson.from.buffer(buf)
checkEquals( typeof(mongo.bson.value(b, "test")), "raw")


# test create bson with list
r <- 1:24
dim(r) <- c(3,2,4)
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append.int(buf, "test", r)
b <- mongo.bson.from.buffer(buf)
checkEquals( mongo.bson.value(b, "test"), r)


# test create bson from data frame
age <- c(5, 8)
height <- c(35, 47)
d <- data.frame(age=age,height=height)
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append.object(buf, "table", d)
b <- mongo.bson.from.buffer(buf)
out <- mongo.bson.value(b, "table")
# checkIdentical( out, d )
# ToDo


# test create bson from list
age=18:29
height=c(76.1,77,78.1,78.2,78.8,79.7,79.9,81.1,81.2,81.8,82.8,83.5)
village=data.frame(age=age,height=height)
unclass(village)
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append.object(buf, "village", village)
b <- mongo.bson.from.buffer(buf)
v <- mongo.bson.value(b, "village")
unclass(v)
#checkEquals(unclass(v), unclass(village))
#ToDo

m <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow=3, ncol=3, dimnames=list(c("X","Y","Z"),c("x","y","z")))
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append.object(buf, "mat", m)
b <- mongo.bson.from.buffer(buf)
v <- mongo.bson.value(b, "mat")
#checkEquals(attributes(m),attributes(v))
# ToDo


buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append(buf, "x", 5L)
scope <- mongo.bson.from.buffer(buf)
codew <- mongo.code.w.scope.create("y = x", scope)
print(codew)
checkEquals( class(codew), c("mongo.code.w.scope", "mongo.code" ))

# create monster object
buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append.string(buf, "name", "Gerald")
mongo.bson.buffer.append.int(buf, "age", 48L)
mongo.bson.buffer.append.bool(buf, "True", TRUE)
mongo.bson.buffer.append.double(buf, "ThreePointFive", 3.5)
mongo.bson.buffer.append.long(buf, "YearSeconds", 365.24219 * 24 * 60 * 60)
mongo.bson.buffer.append.time(buf, "lt", strptime("05-12-2012", "%m-%d-%Y"))
oid <- mongo.oid.from.string("1234567890AB1234567890AB")
mongo.bson.buffer.append.oid(buf, "_id", oid)
id <- mongo.oid.create()
print(mongo.oid.time(id))
print(as.character(id))
mongo.bson.buffer.append(buf, "ID", id)
mongo.bson.buffer.append.null(buf, "Null")
mongo.bson.buffer.start.object(buf, "One_Four")
for (x in 1:4)
  mongo.bson.buffer.append.int(buf, as.character(x), x)
mongo.bson.buffer.finish.object(buf)
code <- mongo.code.create("CoDe")
mongo.bson.buffer.append(buf, "code", code)
mongo.bson.buffer.append(buf, "CodeW", codew)
mongo.bson.buffer.append.symbol(buf, "symbol", "SyMbOl")
mongo.bson.buffer.append.undefined(buf, "undefined1")
undef <- mongo.undefined.create()
mongo.bson.buffer.append(buf, "undefined2", undef)
mongo.bson.buffer.append(buf, "regex",
                         mongo.regex.create("pattern", "options"))
bin <- raw(length=3)
for (i in 0:2)
  bin[i] = as.raw(i * 3)
mongo.bson.buffer.append(buf, "bin", bin)
mongo.bson.buffer.append.time(buf, "Now", Sys.time())
ts <- mongo.timestamp.create(Sys.time() + 24 * 60 * 60, 25L)
mongo.bson.buffer.append.timestamp(buf, "Later", ts)
mongo.bson.buffer.start.object(buf, "data")
mongo.bson.buffer.append(buf, "sub1", 1L)
mongo.bson.buffer.append(buf, "sub2", Sys.time())
mongo.bson.buffer.finish.object(buf)
b <- mongo.bson.from.buffer(buf)
print(b)
checkEquals( class(b), "mongo.bson")


buf <- mongo.bson.buffer.create()
mongo.bson.buffer.append(buf, "name", "Fred")
mongo.bson.buffer.append(buf, "city", "Dayton")
mongo.bson.buffer.append(buf, "age", 21L)
out1 <- mongo.bson.buffer.size(buf)
print(out1)
checkTrue(is.integer(out1))

y <- mongo.bson.from.buffer(buf)
out2 <- mongo.bson.size(y)
print(out2)
checkTrue(is.integer(out2))
checkEquals(out1, out2)


l <- mongo.bson.to.list(y)
checkEquals(l$name, "Fred")
checkEquals(l$city, "Dayton")


buf <- mongo.bson.buffer.create()
l <- list(fruit = "apple", hasSeeds = TRUE)
mongo.bson.buffer.append.list(buf, "item", l)
b <- mongo.bson.from.buffer(buf)
print(b)
checkEquals(mongo.bson.value(b, "item.fruit"), "apple")
checkEquals(mongo.bson.value(b, "item.hasSeeds"), TRUE)


buf <- mongo.bson.buffer.create()
undef <- mongo.undefined.create()
mongo.bson.buffer.append(buf, "Undef", undef)
l <- list(u1 = undef, One = 1)
mongo.bson.buffer.append.list(buf, "listWundef", l)
b <- mongo.bson.from.buffer(buf)
print(b)
checkEquals( class(b), "mongo.bson")
checkEquals( class( mongo.bson.value(b, "Undef") ), "mongo.undefined")
checkEquals( class( mongo.bson.value(b, "listWundef.u1") ), "mongo.undefined")

out <- mongo.bson.to.list(b)
print(out)
checkEquals( class(out$Undef), "mongo.undefined")


# check mongo.bson.to.list
json <- '{"state":"AL"}'
bson <- mongo.bson.from.JSON(json)
list <- mongo.bson.to.list(bson)
checkTrue(is.list(list))
json2 <- toJSON(list)
# checkEquals(json, json2)
