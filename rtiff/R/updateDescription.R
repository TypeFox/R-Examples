"updateDescription" <-
function(fn, description) {
    .Call("updateTTag", fn, description, PACKAGE="rtiff")
    return()
}
