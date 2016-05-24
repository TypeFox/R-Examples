## ----eval = FALSE--------------------------------------------------------
#  modify_lang <- function(x, f, ...) {
#    recurse <- function(y) {
#      # if (!is.null(names(y))) names(y) <- f2(names(y))
#      lapply(y, modify_lang, f = f, ...)
#    }
#  
#    if (is.atomic(x) || is.name(x)) {
#      # Leaf
#      f(x, ...)
#    } else if (is.call(x)) {
#      as.call(recurse(x))
#    } else if (is.function(x)) {
#      formals(x) <- modify_lang(formals(x), f, ...)
#      body(x) <- modify_lang(body(x), f, ...)
#      x
#    } else if (is.pairlist(x)) {
#      # Formal argument lists (when creating functions)
#      as.pairlist(recurse(x))
#    } else if (is.expression(x)) {
#      # shouldn't occur inside tree, but might be useful top-level
#      as.expression(recurse(x))
#    } else if (is.list(x)) {
#      # shouldn't occur inside tree, but might be useful top-level
#      recurse(x)
#    } else {
#      stop("Unknown language class: ", paste(class(x), collapse = "/"),
#        call. = FALSE)
#    }
#  }

## ------------------------------------------------------------------------
identical({ 1 + 2; 3 + 4 }, `{`(1 + 2, 3 + 4))

## ---- eval = FALSE-------------------------------------------------------
#  `{`(count(), as.call(recurse(x)))

## ------------------------------------------------------------------------
f1 <- function() 1

f1 <- function() 2
f1() == 2

## ------------------------------------------------------------------------
env <- new.env()
f1 <- function() 1
env$f2 <- function() f1() + 1

env$f1 <- function() 2

env$f2() == 3

## ----eval = FALSE--------------------------------------------------------
#  replacements_S4 <- function(env) {
#    generics <- getGenerics(env)
#  
#    unlist(recursive = FALSE,
#      Map(generics@.Data, generics@package, USE.NAMES = FALSE,
#        f = function(name, package) {
#        what <- methodsPackageMetaName("T", paste(name, package, sep = ":"))
#  
#        table <- get(what, envir = env)
#  
#        lapply(ls(table, all.names = TRUE), replacement, env = table)
#      })
#    )
#  }

