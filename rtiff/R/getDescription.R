"getDescription" <-
function(fn) {
    ret <- .Call("getTiffDescription", fn, PACKAGE="rtiff")
    ret
}
