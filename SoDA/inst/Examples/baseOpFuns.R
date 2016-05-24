maybeOps <- objects("package:base", all.names=TRUE)
nameRegexp <- "^[.[:alpha:]][._[:alnum:]]*$"
maybeOps <- maybeOps[-grep(nameRegexp, maybeOps)]
nameGetsRegexp <-  "^[.[:alpha:]][._[:alnum:]]*<-$"
maybeOps <- maybeOps[-grep(nameGetsRegexp, maybeOps)]
S3MethodRegexp <- "[.][[:alpha:]][._[:alnum:]]*$"
maybeOps <- maybeOps[-grep(S3MethodRegexp, maybeOps)]
maybeOps
