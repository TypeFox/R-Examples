message("TESTING: trim()...")

library("R.oo")

x <- character(0L)
y <- trim(x)
print(y)
stopifnot(identical(y, x))

x <- "Hello world!"
y <- trim(x)
print(y)
stopifnot(identical(y, x))

x <- " \tHello world!\n "
y <- trim(x)
print(y)
stopifnot(identical(y, "Hello world!"))

x <- c(" \tHello",  "world!")
y <- trim(x)
print(y)
stopifnot(identical(y, c("Hello", "world!")))


message("TESTING: trim()...DONE")
