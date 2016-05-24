### R code from vignette source 'RProtoBuf-quickref.Rnw'

###################################################
### code chunk number 1: RProtoBuf-quickref.Rnw:26-29
###################################################
options(width= 50)
library("RProtoBuf")
rpb.version <- packageDescription("RProtoBuf")$Version


###################################################
### code chunk number 2: RProtoBuf-quickref.Rnw:42-45
###################################################
ab.proto <- system.file( "proto", "addressbook.proto",
	package = "RProtoBuf" )
writeLines( readLines( ab.proto ) )


###################################################
### code chunk number 3: RProtoBuf-quickref.Rnw:50-53 (eval = FALSE)
###################################################
## readProtoFiles( "somefile.proto" )
## readProtoFiles( dir = somedir )
## readProtoFiles( package = AnRPackage )


###################################################
### code chunk number 4: RProtoBuf-quickref.Rnw:58-61
###################################################
message <- new( tutorial.Person, id = 0,
	name = "Romain Francois",
	email = "francoisromain@free.fr" )


###################################################
### code chunk number 5: RProtoBuf-quickref.Rnw:66-77
###################################################
# to a file
tf1 <- tempfile()
message$serialize( tf1 )

# to a binary connection
tf2 <- tempfile(); con <- file( tf2, open = "wb" )
message$serialize( con )
close(con)

# retrieve the payload
message$serialize( NULL )


###################################################
### code chunk number 6: RProtoBuf-quickref.Rnw:82-87
###################################################
# from a file
tutorial.Person$read( tf1 )
# from a connection
con <- file( tf2, open = "rb" )
tutorial.Person$read( con )


###################################################
### code chunk number 7: RProtoBuf-quickref.Rnw:88-89
###################################################
close( con )


###################################################
### code chunk number 8: RProtoBuf-quickref.Rnw:94-98
###################################################
email <- message$email
message$id <- 2
message[[ "name" ]] <- "Romain"
id <- message[[ 2 ]] # tag number for 'id'


###################################################
### code chunk number 9: RProtoBuf-quickref.Rnw:130-131
###################################################
writeLines( message$toString() )


