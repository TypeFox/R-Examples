proto <- function (. = parent.env(envir), expr = {}, envir = 
		new.env(parent = parent.frame()), ..., funEnvir = envir) {
    parent.env(envir) <- .
    as.proto.environment(envir)  # must do this before eval(...)
    # moved eval after for so that ... always done first
    # eval(substitute(eval(quote({ expr }))), envir)
    dots <- list(...); names <- names(dots)
    for (i in seq(length = length(dots))) { 
      assign(names[i], dots[[i]], envir = envir)
      if (!identical(funEnvir, FALSE) && is.function(dots[[i]])) 
        environment(envir[[names[i]]]) <- funEnvir
    }
    eval(substitute(eval(quote({ expr }))), envir)
    if (length(dots)) as.proto.environment(envir) else envir
}

as.proto <- function(x, ...) UseMethod("as.proto")
as.proto.environment <- function(x, ...) {
	assign(".that", x, envir = x)
	assign(".super", parent.env(x), envir = x)
	structure(x, class = c("proto", "environment"))
}
as.proto.proto <- function(x, ...) x
as.proto.list <- function(x, envir, parent, all.names = FALSE, ..., 
	funEnvir = envir, SELECT = function(x) TRUE) {
       if (missing(envir)) {
		if (missing(parent)) parent <- parent.frame()
		envir <- if (is.proto(parent)) 
			parent$proto(...) 
		else 
			proto(parent, ...)
       }
       for(s in names(x))
          if (SELECT(x[[s]])) {
             assign(s, x[[s]], envir = envir)
             if (is.function(x[[s]]) && !identical(funEnvir, FALSE)) 
		environment(envir[[s]]) <- funEnvir
          }
       if (!missing(parent)) parent.env(envir) <- parent
       as.proto.environment(envir)  # force refresh of .that and .super
}

"$<-.proto" <- function(this,s,value) { 
        if (s == ".super") parent.env(this) <- value
	if (is.function(value))  environment(value) <- this
	this[[as.character(substitute(s))]] <- value
	this
}

is.proto <- function(x) inherits(x, "proto")
isnot.function <- function(x) ! is.function(x)

"$.proto" <-
function (this, x, args) {
   inh <- substr(x, 1, 2) != ".."
   p <- parent.frame()
   res <- get(x, envir = this, inherits = inh)
   is.function <- is.function(res)
   is.that <- match(deparse(substitute(this)), c(".that", ".super"),
       nomatch = 0)
   if (is.function && !is.that) {
       res <- function(...) get(x, envir = this, inherits = inh)(this, ...)
       class(res) <- c("instantiatedProtoMethod", "function")
       attr(res, "this") <- this
       if (!missing(args)) res <- do.call(res, args, envir = p)
   }
   res
}


# modified from Tom Short's original
print.instantiatedProtoMethod <- function(x, ...) {
  # cat("proto method call: ")
  # print(unclass(x))
  cat("proto method (instantiated with ", name.proto(attr(x, "this")), 
    "): ", sep = "")
  print(eval(body(x)[[1]], envir = environment(x)))
}

# modified from Tom Short's original
str.proto <- function(object, max.level = 1, nest.lev = 0, indent.str = 
   paste(rep.int(" ", max(0, nest.lev + 1)), collapse = ".."), ...) {
 cat("proto", name.proto(object), "\n")
 Lines <- capture.output(str(as.list(object), max.level = max.level, 
    nest.lev = nest.lev, ...))[-1]
 for(s in Lines) cat(s, "\n")
 if (is.proto(parent.env(object))) {
   cat(indent.str, "parent: ", sep = "")
   str(parent.env(object), nest.lev = nest.lev + 1, ...)
 }
}

