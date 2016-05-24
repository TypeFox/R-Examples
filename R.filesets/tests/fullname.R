library("R.filesets")

message("*** fullname and friends")

name <- "Name"
tags <- c("tag1", "tag2", "tag3")
fullname <- paste(c(name, tags), collapse=",")


message("- fullname()")
fn <- fullname(fullname)
print(fn)
stopifnot(fn == fullname)

fn <- fullname(name, tags=tags)
print(fn)
stopifnot(fn == fullname)

fn <- fullname(name, tags=NULL, tags)
print(fn)
stopifnot(fn == fullname)

fn <- fullname(c(name, tags))
print(fn)
stopifnot(fn == fullname)

parts <- fullname(name, tags=tags, collapse=FALSE)
print(parts)
stopifnot(all(parts == c(name, tags)))

parts <- fullname(fullname, collapse=FALSE)
print(parts)
stopifnot(all(parts == c(name, tags)))

## Argument 'name' is missing
res <- try(fullname(tags=c(name, tags)), silent=TRUE)
stopifnot(inherits(res, "try-error"))


message("- name()")
n <- name(fullname)
print(n)
stopifnot(n == name)

n <- name(parts)
print(n)
stopifnot(n == name)

n <- name(name, tags)
print(n)
stopifnot(n == name)

n <- name(paste(c(name, tags[-1]), collapse=","))
print(n)
stopifnot(n == name)


message("- tags()")
ts <- tags(fullname)
print(ts)
stopifnot(all(ts == tags))

ts <- tags(fullname, collapse=TRUE)
print(ts)
stopifnot(all(ts == paste(tags, collapse=",")))

ts <- tags(parts)
print(ts)
stopifnot(all(ts == tags))

ts <- tags(name, tags)
print(ts)
stopifnot(all(ts == tags))

ts <- tags(paste(c(name, tags[-1]), collapse=","))
print(ts)
stopifnot(all(ts == tags[-1]))


message("- dropTags()")
fn <- dropTags(fullname)
print(fn)
stopifnot(all(fn == parts))

fn <- dropTags(parts)
print(fn)
stopifnot(all(fn == parts))

fn <- dropTags(name, tags)
print(fn)
stopifnot(all(fn == parts))

fn <- dropTags(paste(c(name, tags[-1]), collapse=","))
print(fn)
stopifnot(all(fn == c(name, tags[-1])))

fn <- dropTags(fullname, drop=NULL)
print(fn)
stopifnot(all(fn == parts))

fn <- dropTags(fullname, drop="foo")
print(fn)
stopifnot(all(fn == parts))

fn <- dropTags(fullname, drop=tags)
print(fn)
stopifnot(all(fn == parts[1]))

fn <- dropTags(fullname, drop=tags[1])
print(fn)
stopifnot(all(fn == parts[-2]))

fn <- dropTags(fullname, drop=name)
print(fn)
stopifnot(all(fn == parts))


message("- Arguments$getTags()")
library("R.utils")
ts <- Arguments$getTags(NULL)
stopifnot(length(ts) == 0, is.null(ts))

ts <- Arguments$getTags(character(0L))
stopifnot(length(ts) == 0, is.null(ts))

ts <- Arguments$getTags(NA_character_)
stopifnot(length(ts) == 0, is.null(ts))

ts <- Arguments$getTags(tags)
print(ts)
stopifnot(all.equal(ts, paste(tags, collapse=",")))

ts <- Arguments$getTags(tags, collapse=NULL)
print(ts)
stopifnot(all.equal(ts, tags))


message("*** fullname and friends ... DONE")

