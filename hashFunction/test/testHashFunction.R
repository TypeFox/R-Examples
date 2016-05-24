library(hashFunction)
x <- 1
print( spooky.32(as.integer(x)))
print( spooky.32(as.character(x)))

print( cityhash.64(as.integer(x)))
print( cityhash.64(as.character(x)))

print( murmur3.32(as.integer(x)))
print( murmur3.32(as.character(x)))

x <- "hello"
print( murmur3.32(as.character(x))) ## should equal to 613153351

