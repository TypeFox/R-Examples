
library("arulesSequences")

### ceeboo 2008, 2016

data(zaki)

z <- as(zaki, "timedsequences")
as(as(z, "sequences"), "data.frame")

s <- similarity(z)
s
all(s == similarity(z, z))

as(s, "dist")

similarity(z, strict = TRUE) - s

similarity(z, method = "dice") - s
similarity(z, method = "cosine") - s

similarity(z, method = "subset", strict = TRUE)
similarity(z, method = "subset")

is.subset(z)
is.subset(z, proper = TRUE)

is.subset(z[3], z)
is.subset(z[3], z, proper = TRUE)

is.superset(z)

is.superset(z[1], z)

###
