str.set <-
function(object, ...)
{
    str(unclass(object), ...)
    writeLines(" - attr(*, \"class\")= chr \"set\"")
}

str.gset <-
function(object, ...)
{
    str(unclass(object), ...)
    writeLines(" - attr(*, \"class\")= chr \"gset\"")
}

str.cset <-
function(object, ...)
{
    str(unclass(object), ...)
    writeLines(" - attr(*, \"class\")= chr \"cset\"")
}
