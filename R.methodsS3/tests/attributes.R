library("R.methodsS3")

message("TESTING: attributes()...")

export <- R.methodsS3:::export
`export<-` <- R.methodsS3:::`export<-`
noexport <- R.methodsS3:::noexport
`S3class<-` <- R.methodsS3:::`S3class<-`


foo <- function() NULL
str(foo)

foo <- export(foo)
str(foo)

export(foo) <- TRUE
str(foo)

foo <- noexport(foo)
str(foo)

foo.Bar <- function(...) NULL
S3class(foo.Bar) <- "Bar"
str(foo)

message("TESTING: attributes()...DONE")
