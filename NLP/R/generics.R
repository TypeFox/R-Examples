content <-
function(x)
    UseMethod("content", x)

`content<-` <-
function(x, value)
    UseMethod("content<-", x)

meta <-
function(x, tag = NULL, ...)
    UseMethod("meta", x)

`meta<-` <-
function(x, tag = NULL, ..., value)
    UseMethod("meta<-", x)

