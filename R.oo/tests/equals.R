message("TESTING: equals()...")

library("R.oo")

print(equals(1,1))
print(equals(1,2))
a <- 1:100
b <- 1:100
print(equals(a,b))

obj <- Object()
print(equals(obj, 1))
print(equals(obj, obj))


message("TESTING: equals()...DONE")
