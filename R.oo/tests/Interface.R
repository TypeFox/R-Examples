library("R.oo")

message("TESTING: Interface...")

int <- Interface()
print(int)

setConstructorS3("MyInterface", function(...) {
  extend(Interface(), "MyInterface")
})
int <- MyInterface()
print(int)

setConstructorS3("MyClass", function(...) {
  extend(Object(), c("MyClass", uses(MyInterface(), Interface())))
})
obj <- MyClass()
print(obj)


## For some reason uses(<string>) doesn't play well with covr
if (!"covr" %in% loadedNamespaces()) {
  setConstructorS3("MyClass", function(...) {
    extend(Object(), c("MyClass", uses("MyInterface", "Interface")))
  })
  obj <- MyClass()
  print(obj)
}

message("TESTING: Interface...DONE")
