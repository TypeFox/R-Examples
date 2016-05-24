message("TESTING: local()...")

library("R.oo")

setConstructorS3("Foo", function() {
  extend(Object(), "Foo")
})

setMethodS3("finalize", "Foo", function(this, ...) {
  cat("Finalized Foo\n")
})

x <- Foo()
print(x)

# Trigger finalizer via garbage collection
rm(list="x")
gc()



local({

setConstructorS3("Bar", function() {
  extend(Object(), "Bar")
})

setMethodS3("finalize", "Bar", function(this, ...) {
  cat("Finalized Bar\n")
})

x <- Bar()
print(x)

# Trigger finalizer via garbage collection
rm(list="x")
gc()

})

message("TESTING: local()...DONE")
