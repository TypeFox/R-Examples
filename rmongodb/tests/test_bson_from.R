library(rmongodb)
library(jsonlite)
library(RUnit)

# 10 tests
# 29.01.2014

out <- mongo.bson.from.JSON('{"name" : "Peter"}')
checkEquals(class(out), "mongo.bson")
checkEquals(mongo.bson.value(out, "name"), "Peter")

out <- mongo.bson.from.JSON('{"name" : {"firstname":"Peter"}, "age":12, "visa":[123,321] }')
checkEquals(class(out), "mongo.bson")
checkEquals(mongo.bson.value(out, "name.firstname"), "Peter")
checkEquals(mongo.bson.value(out, "age"), 12)
checkEquals(as.vector( mongo.bson.value(out, "visa") ), c(123,321))

checkException( mongo.bson.from.JSON( "{'name': 'Peter'}", "Not a valid JSON content: {'name': 'Peter'}"))

json <- '{"a":{"b":[1,{"a":3},3]}}'
cat(prettify(json))
validate(json)
out <- mongo.bson.from.JSON( json )
out2 <- mongo.bson.to.list(out)
checkEquals(class(out), "mongo.bson")
checkEquals(class(out2), "list")
checkTrue( is.vector(out2$a$b) )

# check bson.from.df
bson_data <- mongo.bson.from.df(cars)
res <- sapply(bson_data, mongo.bson.to.Robject)
checkEquals(cars, as.data.frame(t(res)))

# Dmitriy Selivanov <selivanov.dmitriy@gmail.com>
# 2014-10-02
# check bson with raw data
a_lst <- list(x="foo", y="bar")
a_bson_raw <- mongo.bson.from.list(lapply(X = a_lst, FUN = charToRaw))
a_simplify_F <- rawToChar(mongo.bson.to.list(b = a_bson_raw, simplify = F)$y)
checkEquals(a_lst$y, a_simplify_F)
a_simplify_T <- rawToChar(mongo.bson.to.list(b = a_bson_raw, simplify = T)$y)
checkEquals(a_lst$y, a_simplify_T)

# also check arrays
lb <- list(x=list(charToRaw("foo"), charToRaw("bar")))
b <- mongo.bson.from.list(lb)
checkEquals(unlist(mongo.bson.to.list(b = b, simplify = T), use.names = F), unlist(lb, use.names = F))
checkEquals(unlist(mongo.bson.to.list(b = b, simplify = F), use.names = T), unlist(lb, use.names = T))

