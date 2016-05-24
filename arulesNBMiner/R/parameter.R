
## helper
.list2object <-  function(from, to) {
    if (!length(from)) return(new(to))
    n = names(from)
    s = slotNames(to)
    p = pmatch(n, s)
    if(any(is.na(p)))
    stop(paste("\nInvalid slot name(s) for class",
            to, ":", paste(n[is.na(p)], collapse=" ")))
    names(from) = s[p]
    do.call("new", c(from, Class=to))
}


## coercion
setAs("NULL", "NBMinerControl", function(from, to) new(to) )
setAs("list", "NBMinerControl", function(from, to) .list2object(from, to))

setAs("NULL", "NBMinerParameter", function(from, to) new(to))
setAs("list", "NBMinerParameter", function(from, to) .list2object(from, to))

## show
setMethod("show", signature(object = "NBMinerControl"),
    function(object) {
        print(data.frame(sapply(slotNames(object),
                    function(x) slot(object, x), simplify = FALSE),
                row.names = ""))
        invisible(object)
    })


setMethod("show", signature(object = "NBMinerParameter"),
    function(object) {
        print(data.frame(sapply(slotNames(object),
                    function(x) slot(object, x),
                    simplify = FALSE), row.names = ""))

        invisible(object)
    })





