### R code from vignette source 'DSL.Rnw'

###################################################
### code chunk number 1: init
###################################################
options(width = 60)
require("DSL")


###################################################
### code chunk number 2: ds_create
###################################################
ds <- DStorage( type = "LFS", base_dir = tempdir(),
                chunksize = 10 * 1024^2 )


###################################################
### code chunk number 3: ds_print_summary
###################################################
ds
summary(ds)


###################################################
### code chunk number 4: dl_create
###################################################
dl <- DList( letters = letters, numbers = 0:9  )
l <- as.list( letters )
names(l) <- LETTERS
dl2 <- as.DList(l)
identical( as.list(dl2), l )
dl3 <- as.DList( system.file("examples", package = "DSL") )


###################################################
### code chunk number 5: dl_create2
###################################################
dl <- DList( letters = letters, numbers = 0:9, DStorage = ds )


###################################################
### code chunk number 6: dl_methods
###################################################
#dl
summary(dl)
names( dl2 )
length(dl3)
dl3[[1]]


###################################################
### code chunk number 7: dl_mapreduce
###################################################
dl <- DList( line1 = "This is the first line.",
             line2 = "Now, the second line." )
res <- DLapply( dl, function(x) unlist(strsplit(x, " ")) )
as.list( res )

foo <- function( keypair )
    list( key = paste("next_", keypair$key, sep = ""), value =
         gsub("first", "mapped", keypair$value) )

dlm <- DMap( x = dl, MAP = foo)
## retrieve keys
unlist(DGather(dlm, keys = TRUE, names = FALSE))
## retrieve values
as.list( dlm )


###################################################
### code chunk number 8: dl_stor_replace (eval = FALSE)
###################################################
## l <- list( line1 = "This is the first line.",
##            line2 = "Now, the second line." )
## dl <- as.DList( l )
## DL_storage(dl)
## ds <- DStorage("HDFS", tempdir())
## DL_storage(dl) <- ds
## as.list(dl)


###################################################
### code chunk number 9: ex1_files
###################################################
## simple wordcount based on two files:
dir(system.file("examples", package = "DSL"))


###################################################
### code chunk number 10: ex1_stor
###################################################
## first force 1 chunk per file (set max chunk size to 1 byte):
ds <- DStorage("LFS", tempdir(), chunksize = 1L)
## make "DList", i.e., read file contents and store in chunks
dl <- as.DList( system.file("examples", package = "DSL"),
                DStorage = ds )


###################################################
### code chunk number 11: ex1_read
###################################################
## read files
dl <- DMap(dl, function( keypair ){
    list( key = keypair$key,
          value = tryCatch(readLines(keypair$value),
                           error = function(x) NA) )
})


###################################################
### code chunk number 12: ex1_map
###################################################
## split into terms
splitwords <- function( keypair ){
    keys <- unlist(strsplit(keypair$value, " "))
    mapply( function(key, value) list( key = key, value = value),
            keys, rep(1L, length(keys)),
            SIMPLIFY = FALSE, USE.NAMES = FALSE )
}
res <- DMap( dl, splitwords )
as.list(res)


###################################################
### code chunk number 13: ex1_reduce
###################################################
## now aggregate by term
res <- DReduce( res, sum )
as.list( res )


