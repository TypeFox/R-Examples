library("R.oo")

message("TESTING: hashCode() ...")

y <- hashCode(character(0L))
print(y)

y <- hashCode("")
print(y)

y <- hashCode(base::letters)
print(y)

y <- hashCode(integer(0L))
print(y)

y <- hashCode(1:10)
print(y)

y <- hashCode(double(0L))
print(y)

y <- hashCode(1:10+0.1)
print(y)

y <- hashCode(list(0L))
print(y)

y <- hashCode(as.list(1:10))
print(y)



message("TESTING: hashCode() ... DONE")
