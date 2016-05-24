### R code from vignette source 'RProtoBuf-intro.Rnw'

###################################################
### code chunk number 1: RProtoBuf-intro.Rnw:26-30
###################################################
library("RProtoBuf")
options("width"=65)
rpb.version <- packageDescription("RProtoBuf")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: RProtoBuf-intro.Rnw:120-121
###################################################
args( readProtoFiles )


###################################################
### code chunk number 3: RProtoBuf-intro.Rnw:127-129
###################################################
proto.dir <- system.file( "proto", package = "RProtoBuf" )
proto.file <- file.path( proto.dir, "addressbook.proto" )


###################################################
### code chunk number 4: RProtoBuf-intro.Rnw:130-131 (eval = FALSE)
###################################################
## readProtoFiles( proto.file )


###################################################
### code chunk number 5: RProtoBuf-intro.Rnw:138-139
###################################################
dir( proto.dir, pattern = "\\.proto$", full.names = TRUE )


###################################################
### code chunk number 6: RProtoBuf-intro.Rnw:140-141 (eval = FALSE)
###################################################
## readProtoFiles( dir = proto.dir )


###################################################
### code chunk number 7: RProtoBuf-intro.Rnw:151-152 (eval = FALSE)
###################################################
## readProtoFiles( package = "RProtoBuf" )


###################################################
### code chunk number 8: RProtoBuf-intro.Rnw:160-161
###################################################
ls( "RProtoBuf:DescriptorPool" )


###################################################
### code chunk number 9: RProtoBuf-intro.Rnw:173-174
###################################################
p <- new( tutorial.Person, name = "Romain", id = 1 )


###################################################
### code chunk number 10: RProtoBuf-intro.Rnw:183-186
###################################################
p$name
p$id
p$email <- "francoisromain@free.fr"


###################################################
### code chunk number 11: RProtoBuf-intro.Rnw:195-198
###################################################
p[["name"]] <- "Romain Francois"
p[[ 2 ]] <- 3
p[[ "email" ]]


###################################################
### code chunk number 12: RProtoBuf-intro.Rnw:213-214
###################################################
p


###################################################
### code chunk number 13: RProtoBuf-intro.Rnw:221-222
###################################################
writeLines( as.character( p ) )


###################################################
### code chunk number 14: RProtoBuf-intro.Rnw:233-234
###################################################
serialize( p, NULL )


###################################################
### code chunk number 15: RProtoBuf-intro.Rnw:239-243
###################################################
tf1 <- tempfile()
tf1
serialize( p, tf1 )
readBin( tf1, raw(0), 500 )


###################################################
### code chunk number 16: RProtoBuf-intro.Rnw:248-253
###################################################
tf2 <- tempfile()
con <- file( tf2, open = "wb" )
serialize( p, con )
close( con )
readBin( tf2, raw(0), 500 )


###################################################
### code chunk number 17: RProtoBuf-intro.Rnw:259-265
###################################################
# serialize to a file
p$serialize( tf1 )
# serialize to a binary connection
con <- file( tf2, open = "wb" )
p$serialize( con )
close( con )


###################################################
### code chunk number 18: RProtoBuf-intro.Rnw:275-276
###################################################
args( read )


###################################################
### code chunk number 19: RProtoBuf-intro.Rnw:285-287
###################################################
message <- read( tutorial.Person, tf1 )
writeLines( as.character( message ) )


###################################################
### code chunk number 20: RProtoBuf-intro.Rnw:293-297
###################################################
con <- file( tf2, open = "rb" )
message <- read( tutorial.Person, con )
close( con )
writeLines( as.character( message ) )


###################################################
### code chunk number 21: RProtoBuf-intro.Rnw:302-305
###################################################
# reading the raw vector payload of the message
payload <- readBin( tf1, raw(0), 5000 )
message <- read( tutorial.Person, payload )


###################################################
### code chunk number 22: RProtoBuf-intro.Rnw:312-320
###################################################
# reading from a file
message <- tutorial.Person$read( tf1 )
# reading from a binary connection
con <- file( tf2, open = "rb" )
message <- tutorial.Person$read( con )
close( con )
# read from the payload
message <- tutorial.Person$read( payload )


###################################################
### code chunk number 23: RProtoBuf-intro.Rnw:331-332
###################################################
str( p )


###################################################
### code chunk number 24: RProtoBuf-intro.Rnw:416-427
###################################################
message <- new( tutorial.Person,
	name = "foo", email = "foo@bar.com", id = 2,
	phone = list(
		new( tutorial.Person.PhoneNumber, number = "+33(0)...", type = "HOME" ),
		new( tutorial.Person.PhoneNumber, number = "+33(0)###", type = "MOBILE" )
	) )
message$name
message$email
message[[ "phone" ]]
# using the tag number
message[[ 2 ]] # id


###################################################
### code chunk number 25: RProtoBuf-intro.Rnw:489-495
###################################################
message <- new( tutorial.Person,
	name = "foo", id = 2 )
message$email <- "foo@bar.com"
message[[ "id" ]] <- 2
message[[ 1 ]] <- "foobar"
writeLines( message$as.character() )


###################################################
### code chunk number 26: RProtoBuf-intro.Rnw:544-548
###################################################
message <- new( tutorial.Person, name = "foo" )
message$has( "name" )
message$has( "id" )
message$has( "phone" )


###################################################
### code chunk number 27: RProtoBuf-intro.Rnw:558-563
###################################################
m1 <- new( tutorial.Person, name = "foo" )
m2 <- m1$clone( )
m2$email <- "foo@bar.com"
writeLines( as.character( m1 ) )
writeLines( as.character( m2 ) )


###################################################
### code chunk number 28: RProtoBuf-intro.Rnw:574-578
###################################################
message <- new( tutorial.Person, name = "foo" )
message$isInitialized()
message$id <- 2
message$isInitialized()


###################################################
### code chunk number 29: RProtoBuf-intro.Rnw:587-597
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
tf1 <- tempfile( )
tf1
message$serialize( tf1 )

tf2 <- tempfile( )
tf2
con <- file( tf2, open = "wb" )
message$serialize( con )
close( con )


###################################################
### code chunk number 30: RProtoBuf-intro.Rnw:604-606
###################################################
readBin( tf1, raw(0), 500 )
readBin( tf2, raw(0), 500 )


###################################################
### code chunk number 31: RProtoBuf-intro.Rnw:612-613
###################################################
message$serialize(NULL)


###################################################
### code chunk number 32: RProtoBuf-intro.Rnw:622-630
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
writeLines( as.character( message ) )
message$clear()
writeLines( as.character( message ) )

message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
message$clear( "id" )
writeLines( as.character( message ) )


###################################################
### code chunk number 33: RProtoBuf-intro.Rnw:643-650
###################################################
message <- new( tutorial.Person, name = "foo",
	phone = list(
		new( tutorial.Person.PhoneNumber, number = "+33(0)...", type = "HOME"  ),
		new( tutorial.Person.PhoneNumber, number = "+33(0)###", type = "MOBILE"  )
		) )
message$size( "phone" )
size( message, "phone" )


###################################################
### code chunk number 34: RProtoBuf-intro.Rnw:664-668
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
message$bytesize()
bytesize( message )
length( message$serialize( NULL ) )


###################################################
### code chunk number 35: RProtoBuf-intro.Rnw:677-689
###################################################
message <- new( tutorial.Person, name = "foo",
	phone = list(
		new( tutorial.Person.PhoneNumber, number = "+33(0)...", type = "HOME"  ),
		new( tutorial.Person.PhoneNumber, number = "+33(0)###", type = "MOBILE"  )
		) )
message$swap( "phone", 1, 2 )
writeLines( as.character( message$phone[[1]] ) )
writeLines( as.character( message$phone[[2]] ) )

swap( message, "phone", 1, 2 )
writeLines( as.character( message$phone[[1]] ) )
writeLines( as.character( message$phone[[2]] ) )


###################################################
### code chunk number 36: RProtoBuf-intro.Rnw:699-708
###################################################
message <- new( tutorial.Person, name = "foo",
	phone = list(
		new( tutorial.Person.PhoneNumber, number = "+33(0)...", type = "HOME"  ),
		new( tutorial.Person.PhoneNumber, number = "+33(0)###", type = "MOBILE"  )
		) )
number <- new( tutorial.Person.PhoneNumber,
		number = "+33(0)---", type = "WORK"  )
message$set( "phone", 1, number )
writeLines( as.character( message ) )


###################################################
### code chunk number 37: RProtoBuf-intro.Rnw:717-723
###################################################
message <- new( tutorial.Person, name = "foo",
	phone = list(
		new( tutorial.Person.PhoneNumber, number = "+33(0)...", type = "HOME"  ),
		new( tutorial.Person.PhoneNumber, number = "+33(0)###", type = "MOBILE"  )
		) )
message$fetch( "phone", 1 )


###################################################
### code chunk number 38: RProtoBuf-intro.Rnw:732-744
###################################################
if (!exists("protobuf_unittest.TestAllTypes",
            "RProtoBuf:DescriptorPool")) {
    unittest.proto.file <- system.file("unitTests", "data",
                                       "unittest.proto",
                                       package="RProtoBuf")
    readProtoFiles(file=unittest.proto.file)
}

## Test setting a singular extensions.
test <- new(protobuf_unittest.TestAllExtensions)
test$setExtension(protobuf_unittest.optional_int32_extension,
                  as.integer(1))


###################################################
### code chunk number 39: RProtoBuf-intro.Rnw:753-754
###################################################
test$getExtension(protobuf_unittest.optional_int32_extension)


###################################################
### code chunk number 40: RProtoBuf-intro.Rnw:762-767
###################################################
message <- new( tutorial.Person, name = "foo")
phone <- new( tutorial.Person.PhoneNumber,
	number = "+33(0)...", type = "HOME"  )
message$add( "phone", phone )
writeLines( message$toString() )


###################################################
### code chunk number 41: RProtoBuf-intro.Rnw:776-779
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
message$str()
str( message )


###################################################
### code chunk number 42: RProtoBuf-intro.Rnw:787-790
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
writeLines( message$as.character() )
writeLines( as.character( message ) )


###################################################
### code chunk number 43: RProtoBuf-intro.Rnw:798-801
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
writeLines( message$toString() )
writeLines( toString( message ) )


###################################################
### code chunk number 44: RProtoBuf-intro.Rnw:809-811
###################################################
message <- new( tutorial.Person, name = "foo", email = "foo@bar.com", id = 2 )
as.list( message )


###################################################
### code chunk number 45: RProtoBuf-intro.Rnw:825-831
###################################################
message <- new( tutorial.Person )
update( message,
	name = "foo",
	id = 2,
	email = "foo@bar.com" )
writeLines( message$as.character() )


###################################################
### code chunk number 46: RProtoBuf-intro.Rnw:841-844
###################################################
message <- new( tutorial.Person )
message$descriptor()
descriptor( message )


###################################################
### code chunk number 47: RProtoBuf-intro.Rnw:855-858
###################################################
message <- new( tutorial.Person )
message$fileDescriptor()
fileDescriptor( message )


###################################################
### code chunk number 48: RProtoBuf-intro.Rnw:937-947
###################################################
# field descriptor
tutorial.Person$email

# enum descriptor
tutorial.Person$PhoneType

# nested type descriptor
tutorial.Person$PhoneNumber
# same as
tutorial.Person.PhoneNumber


###################################################
### code chunk number 49: RProtoBuf-intro.Rnw:956-958
###################################################
tutorial.Person$new( )
new( tutorial.Person )


###################################################
### code chunk number 50: RProtoBuf-intro.Rnw:964-968
###################################################
tutorial.Person$new( email = "foo@bar.com" )

# same as
update( tutorial.Person$new( ), email = "foo@bar.com" )


###################################################
### code chunk number 51: RProtoBuf-intro.Rnw:977-989
###################################################
# start by serializing a message
message <- new( tutorial.Person.PhoneNumber,
	type = "HOME", number = "+33(0)...." )
tf <- tempfile()
serialize( message, tf )

# now read back the message
m <- tutorial.Person.PhoneNumber$read( tf )
writeLines( as.character( m ) )

m <- read( tutorial.Person.PhoneNumber, tf )
writeLines( as.character( m ) )


###################################################
### code chunk number 52: RProtoBuf-intro.Rnw:999-1004
###################################################
# start by generating the ASCII representation of a message
text <- as.character(new(tutorial.Person, id=1, name="Murray"))
text
# Then read the ascii representation in as a new message object.
msg <- tutorial.Person$readASCII(text)


###################################################
### code chunk number 53: RProtoBuf-intro.Rnw:1019-1023
###################################################
desc <- tutorial.Person
writeLines( desc$toString() )
writeLines( toString( desc ) )
writeLines( as.character(tutorial.Person) )


###################################################
### code chunk number 54: RProtoBuf-intro.Rnw:1032-1033
###################################################
tutorial.Person$as.list()


###################################################
### code chunk number 55: RProtoBuf-intro.Rnw:1042-1043
###################################################
tutorial.Person$asMessage()


###################################################
### code chunk number 56: RProtoBuf-intro.Rnw:1053-1056
###################################################
desc <- tutorial.Person
desc$fileDescriptor()
fileDescriptor( desc )


###################################################
### code chunk number 57: RProtoBuf-intro.Rnw:1065-1069
###################################################
# simple name
tutorial.Person$name()
# name including scope
tutorial.Person$name(full = TRUE)


###################################################
### code chunk number 58: RProtoBuf-intro.Rnw:1078-1080
###################################################
tutorial.Person$containing_type()
tutorial.Person$PhoneNumber$containing_type()


###################################################
### code chunk number 59: RProtoBuf-intro.Rnw:1089-1090
###################################################
tutorial.Person$field_count()


###################################################
### code chunk number 60: RProtoBuf-intro.Rnw:1099-1100
###################################################
tutorial.Person$field(1)


###################################################
### code chunk number 61: RProtoBuf-intro.Rnw:1109-1110
###################################################
tutorial.Person$nested_type_count()


###################################################
### code chunk number 62: RProtoBuf-intro.Rnw:1119-1120
###################################################
tutorial.Person$nested_type(1)


###################################################
### code chunk number 63: RProtoBuf-intro.Rnw:1129-1130
###################################################
tutorial.Person$enum_type_count()


###################################################
### code chunk number 64: RProtoBuf-intro.Rnw:1139-1140
###################################################
tutorial.Person$enum_type(1)


###################################################
### code chunk number 65: RProtoBuf-intro.Rnw:1213-1214
###################################################
writeLines( as.character( tutorial.Person$PhoneNumber ) )


###################################################
### code chunk number 66: RProtoBuf-intro.Rnw:1222-1223
###################################################
writeLines( tutorial.Person.PhoneNumber$toString() )


###################################################
### code chunk number 67: RProtoBuf-intro.Rnw:1232-1234
###################################################
tutorial.Person$id$asMessage()
writeLines(as.character(tutorial.Person$id$asMessage()))


###################################################
### code chunk number 68: RProtoBuf-intro.Rnw:1243-1247
###################################################
# simple name.
name( tutorial.Person$id )
# name including scope.
name( tutorial.Person$id, full=TRUE )


###################################################
### code chunk number 69: RProtoBuf-intro.Rnw:1256-1258
###################################################
fileDescriptor(tutorial.Person$id)
tutorial.Person$id$fileDescriptor()


###################################################
### code chunk number 70: RProtoBuf-intro.Rnw:1267-1269
###################################################
containing_type(tutorial.Person$id)
tutorial.Person$id$containing_type()


###################################################
### code chunk number 71: RProtoBuf-intro.Rnw:1279-1281
###################################################
is_extension( tutorial.Person$id )
tutorial.Person$id$is_extension()


###################################################
### code chunk number 72: RProtoBuf-intro.Rnw:1289-1291
###################################################
number( tutorial.Person$id )
tutorial.Person$id$number()


###################################################
### code chunk number 73: RProtoBuf-intro.Rnw:1300-1302
###################################################
type( tutorial.Person$id )
tutorial.Person$id$type()


###################################################
### code chunk number 74: RProtoBuf-intro.Rnw:1310-1312
###################################################
cpp_type( tutorial.Person$id )
tutorial.Person$id$cpp_type()


###################################################
### code chunk number 75: RProtoBuf-intro.Rnw:1324-1327
###################################################
label( tutorial.Person$id )
label( tutorial.Person$id , TRUE)
tutorial.Person$id$label(TRUE)


###################################################
### code chunk number 76: RProtoBuf-intro.Rnw:1334-1336
###################################################
is_repeated( tutorial.Person$id )
tutorial.Person$id$is_repeated()


###################################################
### code chunk number 77: RProtoBuf-intro.Rnw:1344-1346
###################################################
is_required( tutorial.Person$id )
tutorial.Person$id$is_required()


###################################################
### code chunk number 78: RProtoBuf-intro.Rnw:1354-1356
###################################################
is_optional( tutorial.Person$id )
tutorial.Person$id$is_optional()


###################################################
### code chunk number 79: RProtoBuf-intro.Rnw:1365-1367
###################################################
has_default_value(tutorial.Person$PhoneNumber$type)
has_default_value(tutorial.Person$PhoneNumber$number)


###################################################
### code chunk number 80: RProtoBuf-intro.Rnw:1375-1377
###################################################
default_value( tutorial.Person$PhoneNumber$type )
default_value( tutorial.Person$PhoneNumber$number )


###################################################
### code chunk number 81: RProtoBuf-intro.Rnw:1386-1388
###################################################
message_type(tutorial.Person$phone)
tutorial.Person$phone$message_type()


###################################################
### code chunk number 82: RProtoBuf-intro.Rnw:1395-1396
###################################################
enum_type(tutorial.Person$PhoneNumber$type)


###################################################
### code chunk number 83: RProtoBuf-intro.Rnw:1468-1470
###################################################
tutorial.Person$PhoneType$WORK
name(tutorial.Person$PhoneType$value(number=2))


###################################################
### code chunk number 84: RProtoBuf-intro.Rnw:1479-1480
###################################################
as.list( tutorial.Person$PhoneType )


###################################################
### code chunk number 85: RProtoBuf-intro.Rnw:1488-1489
###################################################
writeLines( as.character( tutorial.Person$PhoneType ) )


###################################################
### code chunk number 86: RProtoBuf-intro.Rnw:1497-1498
###################################################
writeLines( toString( tutorial.Person$PhoneType ) )


###################################################
### code chunk number 87: RProtoBuf-intro.Rnw:1507-1509
###################################################
tutorial.Person$PhoneType$asMessage()
writeLines(as.character(tutorial.Person$PhoneType$asMessage()))


###################################################
### code chunk number 88: RProtoBuf-intro.Rnw:1518-1522
###################################################
# simple name.
name( tutorial.Person$PhoneType )
# name including scope.
name( tutorial.Person$PhoneType, full=TRUE )


###################################################
### code chunk number 89: RProtoBuf-intro.Rnw:1531-1533
###################################################
fileDescriptor(tutorial.Person$PhoneType)
tutorial.Person$PhoneType$fileDescriptor()


###################################################
### code chunk number 90: RProtoBuf-intro.Rnw:1542-1543
###################################################
tutorial.Person$PhoneType$containing_type()


###################################################
### code chunk number 91: RProtoBuf-intro.Rnw:1551-1553
###################################################
length(tutorial.Person$PhoneType)
tutorial.Person$PhoneType$length()


###################################################
### code chunk number 92: RProtoBuf-intro.Rnw:1562-1564
###################################################
tutorial.Person$PhoneType$has("WORK")
tutorial.Person$PhoneType$has("nonexistant")


###################################################
### code chunk number 93: RProtoBuf-intro.Rnw:1573-1575
###################################################
value_count(tutorial.Person$PhoneType)
tutorial.Person$PhoneType$value_count()


###################################################
### code chunk number 94: RProtoBuf-intro.Rnw:1585-1588
###################################################
tutorial.Person$PhoneType$value(1)
tutorial.Person$PhoneType$value(name="HOME")
tutorial.Person$PhoneType$value(number=1)


###################################################
### code chunk number 95: RProtoBuf-intro.Rnw:1647-1648
###################################################
number( tutorial.Person$PhoneType$value(number=2) )


###################################################
### code chunk number 96: RProtoBuf-intro.Rnw:1657-1661
###################################################
# simple name.
name( tutorial.Person$PhoneType$value(number=2) )
# name including scope.
name( tutorial.Person$PhoneType$value(number=2), full=TRUE )


###################################################
### code chunk number 97: RProtoBuf-intro.Rnw:1670-1671
###################################################
enum_type( tutorial.Person$PhoneType$value(number=2) )


###################################################
### code chunk number 98: RProtoBuf-intro.Rnw:1680-1681
###################################################
writeLines( as.character( tutorial.Person$PhoneType$value(number=2) ) )


###################################################
### code chunk number 99: RProtoBuf-intro.Rnw:1690-1691
###################################################
writeLines( toString( tutorial.Person$PhoneType$value(number=2) ) )


###################################################
### code chunk number 100: RProtoBuf-intro.Rnw:1700-1702
###################################################
tutorial.Person$PhoneType$value(number=2)$asMessage()
writeLines(as.character(tutorial.Person$PhoneType$value(number=2)$asMessage()))


###################################################
### code chunk number 101: RProtoBuf-intro.Rnw:1736-1739
###################################################
f <- tutorial.Person$fileDescriptor()
f
f$Person


###################################################
### code chunk number 102: RProtoBuf-intro.Rnw:1771-1772
###################################################
writeLines( as.character(fileDescriptor(tutorial.Person)) )


###################################################
### code chunk number 103: RProtoBuf-intro.Rnw:1780-1781
###################################################
writeLines( fileDescriptor(tutorial.Person)$toString() )


###################################################
### code chunk number 104: RProtoBuf-intro.Rnw:1789-1791
###################################################
asMessage(tutorial.Person$fileDescriptor())
writeLines( as.character(asMessage(tutorial.Person$fileDescriptor())) )


###################################################
### code chunk number 105: RProtoBuf-intro.Rnw:1799-1800
###################################################
as.list( tutorial.Person$fileDescriptor() )


###################################################
### code chunk number 106: RProtoBuf-intro.Rnw:1810-1812
###################################################
name( tutorial.Person$fileDescriptor() )
tutorial.Person$fileDescriptor()$name(TRUE)


###################################################
### code chunk number 107: RProtoBuf-intro.Rnw:1821-1822
###################################################
tutorial.Person$fileDescriptor()$package()


###################################################
### code chunk number 108: RProtoBuf-intro.Rnw:1848-1860
###################################################
# coerce a message type descriptor to a message
# asMessage( tutorial.Person )

# coerce a enum descriptor
asMessage( tutorial.Person.PhoneType )

# coerce a field descriptor
asMessage( tutorial.Person$email )

# coerce a file descriptor
asMessage( fileDescriptor( tutorial.Person ) )



###################################################
### code chunk number 109: RProtoBuf-intro.Rnw:1891-1901
###################################################
message <- new( tutorial.Person, email = "foo@bar.com" )
with( message, {
	# set the id field
	id <- 2

	# set the name field from the email field
	name <- gsub( "[@]", " ", email )

	sprintf( "%d [%s] : %s", id, email, name )
} )


###################################################
### code chunk number 110: RProtoBuf-intro.Rnw:1913-1916
###################################################
m1 <- new( tutorial.Person, email = "foo@bar.com", id = 2 )
m2 <- update( new( tutorial.Person) , email = "foo@bar.com", id = 2 )
identical( m1, m2 )


###################################################
### code chunk number 111: RProtoBuf-intro.Rnw:1921-1923
###################################################
m1 == m2
m1 != m2


###################################################
### code chunk number 112: RProtoBuf-intro.Rnw:1933-1937
###################################################
m1 <- new( tutorial.Person, name = "foobar" )
m2 <- new( tutorial.Person, email = "foo@bar.com" )
m3 <- merge( m1, m2 )
writeLines( as.character( m3 ) )


###################################################
### code chunk number 113: RProtoBuf-intro.Rnw:1946-1952
###################################################
P("tutorial.Person")
new( P("tutorial.Person") )

# but we can do this instead
tutorial.Person
new( tutorial.Person )


###################################################
### code chunk number 114: RProtoBuf-intro.Rnw:1977-1984
###################################################
  extend.proto <- tempfile()
  writeLines(c(
               paste0('import "',
                      name(tutorial.Person$fileDescriptor(), TRUE), '";'),
               "package tutorial;",
               paste0("extend Person {\n  optional string nationality = 100;\n}")),
             extend.proto)


###################################################
### code chunk number 115: RProtoBuf-intro.Rnw:1987-1988
###################################################
writeLines(readLines(extend.proto))


###################################################
### code chunk number 116: RProtoBuf-intro.Rnw:1995-2002
###################################################
library(RProtoBuf)
readProtoFiles(extend.proto)
person <- new(tutorial.Person, id=1, name="Murray")
person
person$setExtension(P("tutorial.nationality"), "USA")
cat(as.character(person))
person$getExtension(P("tutorial.nationality"))


###################################################
### code chunk number 117: RProtoBuf-intro.Rnw:2034-2035
###################################################
2^53 == (2^53 + 1)


###################################################
### code chunk number 118: RProtoBuf-intro.Rnw:2044-2051
###################################################
if (!exists("protobuf_unittest.TestAllTypes",
            "RProtoBuf:DescriptorPool")) {
    unittest.proto.file <- system.file("unitTests", "data",
                                       "unittest.proto",
                                       package="RProtoBuf")
    readProtoFiles(file=unittest.proto.file)
}


###################################################
### code chunk number 119: RProtoBuf-intro.Rnw:2150-2151
###################################################
options("RProtoBuf.int64AsString" = FALSE)


###################################################
### code chunk number 120: RProtoBuf-intro.Rnw:2170-2174
###################################################
test <- new(protobuf_unittest.TestAllTypes)
test$optionalgroup$a <- 3
test$optionalgroup$a
cat(as.character(test))


