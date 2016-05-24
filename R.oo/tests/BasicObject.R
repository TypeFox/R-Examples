message("TESTING: BasicObject...")

library("R.oo")

obj <- BasicObject()
print(obj)

obj <- BasicObject(42L)
print(obj)
stopifnot(obj == 42L)

obj$a <- 1:99
print(obj)

fields <- getFields(obj)
print(fields)

hasA <- hasField(obj, "a")
print(hasA)
stopifnot(hasA)

value <- obj$a
str(value)
stopifnot(identical(value, 1:99))

obj$a <- 1:100
print(obj)

value <- obj[["a"]]
str(value)
stopifnot(identical(value, 1:100))

size <- objectSize(obj)
print(size)

ref <- isReferable(obj)
print(ref)
stopifnot(isTRUE(ref))

time <- getInstantiationTime(obj)
print(time)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Attach and detach
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
obj <- BasicObject(42L)
obj$a <- 1:99

res <- attach(obj)
print(res)
stopifnot(exists("a", mode="integer"))
str(a)

## Object already attached
res <- tryCatch(attach(obj), warning=function(w) w)
stopifnot(inherits(res, "warning"))

res <- detach(obj)
print(res)

## Object already detached
res <- tryCatch(detach(obj), warning=function(w) w)
stopifnot(inherits(res, "warning"))


obj <- BasicObject(list(a=1L, b=2, c=3))

res <- attach(obj)
print(res)
stopifnot(exists("a", mode="integer"))
str(a)

res <- detach(obj)
print(res)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Class
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
obj2 <- newInstance(obj, 43L)
print(obj2)
stopifnot(obj2 == 43L)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Inheritance
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setConstructorS3("MyObject", function(...) {
  extend(BasicObject(), "MyObject", ...)
})
obj <- MyObject(a=1, b=2)
print(obj)
str(obj)
stopifnot(all(c("a", "b") %in% names(attributes(obj))))


setMethodS3("foo", "MyObject", function(static, x=1L, ...) {
  list(x=x, ...)
}, static=TRUE)

res <- MyObject$foo(y=2L)
stopifnot(identical(res$x, 1L))
stopifnot(identical(res$y, 2L))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FIXME: hashCode() return integer(0) whenever
# getInstantiationTime() returns NULL, which is
# now the default behavior of BasicObject
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hash <- hashCode(obj)
print(hash)
## FIXME: Currently returns integer(0)
## stopifnot(length(hash) == 1L)

neq <- equals(obj, 1)
print(neq)
## FIXME: Currently returns NA
## stopifnot(!neq)

eq <- equals(obj, obj)
print(eq)
## FIXME: Currently returns NA
## stopifnot(eq)


message("TESTING: BasicObject...DONE")
