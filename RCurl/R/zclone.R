setGeneric("clone", function(x, ...) standardGeneric("clone"))
setMethod("clone", "ANY", function(x, ...) x)

setMethod("clone", "CURLHandle", function(x, ...) dupCurlHandle(x, ...))
setMethod("clone", "environment",
            function(x, ...) {
               e = new.env(parent = parent.env(x))
               sapply(ls(x, all.names = TRUE),
                      function(id)
                       assign(id, get(id, x), e))
               e
            })
