makeObject = function(.self) {
  force(.self)
  list(
    ls = function(pattern = NULL) {
      Ls(.self, pattern)
    },
    get = function(key, simplify, use.cache) {
      Get(.self, asKeys(.self, key, len = 1L),
        use.cache = asFlag(use.cache, default = .self$use.cache))
    },
    pos = function(n = 1L, use.cache) {
      keys = Ls(.self)
      if (n > length(keys))
        return(NULL)
      Get(.self, keys[n], use.cache = asFlag(use.cache, default = .self$use.cache))
    },
    put = function(..., keys, li = list(), use.cache) {
      Put(.self, ..., keys = keys, li = as.list(li),
        use.cache = asFlag(use.cache, default = .self$use.cache))
    },
    remove = function(keys) {
      Remove(.self, asKeys(.self, keys))
    },
    as.list = function(keys, use.cache) {
      AsList(.self, asKeys(.self, keys, default = Ls(.self)),
        use.cache = asFlag(use.cache, default = .self$use.cache))
    },
    apply = function(FUN, ..., keys, use.cache, simplify = FALSE, use.names = TRUE) {
      Apply(.self, FUN, ..., keys = asKeys(.self, keys, default = Ls(.self)),
        use.cache = asFlag(use.cache, default = .self$use.cache),
        simplify = asFlag(simplify), use.names = asFlag(use.names))
    },
    mapply = function(FUN, ..., keys, use.cache, moreArgs = NULL, simplify = FALSE, use.names = TRUE) {
      Mapply(.self, FUN, ..., keys = asKeys(.self, keys, default = Ls(.self)),
        use.cache = asFlag(use.cache, default = .self$use.cache),
        moreArgs = as.list(moreArgs), simplify = asFlag(simplify), use.names = asFlag(use.names))
    },
    assign = function(keys, envir = parent.frame(), use.cache) {
      Assign(.self, keys = asKeys(.self, keys, default = Ls(.self)), envir = as.environment(envir),
        use.cache = asFlag(use.cache, default = .self$use.cache))
    },
    size = function(keys, unit = "b") {
      match.arg(unit, choices = names(UNITCONVERT))
      Size(.self, asKeys(.self, keys, default = Ls(.self)), unit = unit)
    },
    clear = function(keys) {
      Clear(.self, asKeys(.self, keys, default = Ls(.self)))
    },
    cached = function() {
      Cached(.self)
    },
    info = function() {
      Info(.self)
    }
  )
}
