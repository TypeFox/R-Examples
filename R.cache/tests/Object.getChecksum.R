library("R.cache")

obj <- R.oo::Object()
obj$value <- 42L
print(getChecksum(obj))
