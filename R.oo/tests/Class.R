message("TESTING: Class...")

library("R.oo")

clazz <- Object
print(clazz)
print(as.character(clazz))


message(" - instantiation ...")
obj <- Object()
print(obj)
obj <- newInstance(clazz)
print(obj)
stopifnot(inherits(obj, "Object"))
obj <- clazz[["newInstance"]]()
print(obj)
stopifnot(inherits(obj, "Object"))
message(" - instantiation ... DONE")


message(" - reflection ...")
clazz <- Class$forName("Object")
print(clazz)
print(getKnownSubclasses(clazz))
static <- getStaticInstance(clazz)
print(static)
pkg <- getPackage(clazz)
print(pkg)
stopifnot(pkg == "R.oo")

## Odds and ends
print(isBeingCreated(clazz))

## FIXME: Case should never occur but code allows for it
print(isBeingCreated.Class(obj))
message(" - reflection ... DONE")


message(" - modifiers ...")
print(isAbstract(clazz))
print(isPrivate(clazz))
print(isProtected(clazz))
print(isPublic(clazz))
print(isDeprecated(clazz))
print(isStatic(clazz))  ## TRUE because of Object$load()
message(" - modifiers ... DONE")


message(" - inheritance ...")
setConstructorS3("MyClass", function(...) {
  extend(Object(), "MyClass", ...)
})

obj <- MyClass(a=1, b=2, c=3)
print(obj)
stopifnot(all(c("a", "b", "c") %in% names(obj)))

obj <- newInstance(MyClass, a=1, b=2, c=3)
print(obj)
stopifnot(all(c("a", "b", "c") %in% names(obj)))

message(" - inheritance ... DONE")



message("TESTING: Class...DONE")
