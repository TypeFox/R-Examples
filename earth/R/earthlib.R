# earthlib.R: general purpose routines for the earth package

sos <- function(x) sum(as.vector(x^2)) # sum of squares

# Modify call to refer to the generic e.g. "foo.default" becomes "foo".
# This means that functions like update() call foo() and not foo.default().
# An advantage is that we don't have to export foo.default().

make.call.generic <- function(call, fname)
{
    call[[1]] <- as.name(fname)
    call
}
