library("globals")
opts <- options(warn=1L)

exprs <- list(
  ok   = substitute(function(...) sum(x, ...)),
  warn = substitute(sum(x, ...))
)


message("*** findGlobals() ...")


for (name in names(exprs)) {
  expr <- exprs[[name]]

  message("\n*** codetools::findGlobals():")
  fun <- globals:::asFunction(expr)
  print(fun)
  globals <- codetools::findGlobals(fun)
  print(globals)
  stopifnot(all.equal(globals, c("sum", "x")))

  message("\n*** findGlobals(dotdotdot='ignore'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- findGlobals(expr, dotdotdot="ignore")
  print(globals)
  stopifnot(all.equal(globals, c("sum", "x")))

  message("\n*** findGlobals(dotdotdot='return'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- findGlobals(expr, dotdotdot="return")
  print(globals)
  if (name == "ok") {
    stopifnot(all.equal(globals, c("sum", "x")))
  } else {
    stopifnot(all.equal(globals, c("sum", "x", "...")))
  }

  message("\n*** findGlobals(dotdotdot='warn'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- findGlobals(expr, dotdotdot="warn")
  print(globals)
  if (name == "ok") {
    stopifnot(all.equal(globals, c("sum", "x")))
  } else {
    stopifnot(all.equal(globals, c("sum", "x", "...")))
  }

  message("\n*** findGlobals(dotdotdot='error'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- try(findGlobals(expr, dotdotdot="error"))
  if (name == "ok") {
    stopifnot(all.equal(globals, c("sum", "x")))
  } else {
    stopifnot(inherits(globals, "try-error"))
  }
} # for (name ...)

message("\n*** findGlobals(<exprs>, dotdotdot='return'):")
print(exprs)
globals <- findGlobals(exprs, dotdotdot="return")
print(globals)


message("*** findGlobals() ... DONE")



message("*** globalsOf() ...")

x <- 1:2

for (name in names(exprs)) {
  expr <- exprs[[name]]

  message("\n*** globalsOf(dotdotdot='ignore'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- globalsOf(expr, dotdotdot="ignore")
  print(globals)
  stopifnot(all.equal(names(globals), c("sum", "x")))
  stopifnot(all.equal(globals$sum, base::sum))
  stopifnot(all.equal(globals$x, x))

  message("\n*** globalsOf(dotdotdot='return'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- globalsOf(expr, dotdotdot="return")
  print(globals)
  if (name == "ok") {
    stopifnot(all.equal(names(globals), c("sum", "x")))
  } else {
    stopifnot(all.equal(names(globals), c("sum", "x", "...")))
    stopifnot(!is.list(globals$`...`) && is.na(globals$`...`))
  }
  stopifnot(all.equal(globals$sum, base::sum))
  stopifnot(all.equal(globals$x, x))

  message("\n*** globalsOf(dotdotdot='warn'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- globalsOf(expr, dotdotdot="warn")
  print(globals)
  if (name == "ok") {
    stopifnot(all.equal(names(globals), c("sum", "x")))
  } else {
    stopifnot(all.equal(names(globals), c("sum", "x", "...")))
    stopifnot(!is.list(globals$`...`) && is.na(globals$`...`))
  }
  stopifnot(all.equal(globals$sum, base::sum))
  stopifnot(all.equal(globals$x, x))

  message("\n*** globalsOf(dotdotdot='error'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- try(globalsOf(expr, dotdotdot="error"))
  if (name == "ok") {
    stopifnot(all.equal(names(globals), c("sum", "x")))
    stopifnot(all.equal(globals$sum, base::sum))
    stopifnot(all.equal(globals$x, x))
  } else {
    stopifnot(inherits(globals, "try-error"))
  }
} # for (name ...)

message("\n*** globalsOf(<exprs>, dotdotdot='return'):")
print(exprs)
globals <- globalsOf(exprs, dotdotdot="return")
print(globals)


message("*** globalsOf() ... DONE")


message("*** function(x, ...) globalsOf() ...")

aux <- function(x, ...) {
  args <- list(...)

for (name in names(exprs)) {
  expr <- exprs[[name]]

  message("\n*** globalsOf(dotdotdot='ignore'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- globalsOf(expr, dotdotdot="ignore")
  print(globals)
  stopifnot(all.equal(names(globals), c("sum", "x")))
  stopifnot(all.equal(globals$sum, base::sum))
  stopifnot(all.equal(globals$x, x))

  message("\n*** globalsOf(dotdotdot='return'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- globalsOf(expr, dotdotdot="return")
  print(globals)
  if (name == "ok") {
    stopifnot(all.equal(names(globals), c("sum", "x")))
  } else {
    stopifnot(all.equal(names(globals), c("sum", "x", "...")))
    stopifnot(all.equal(globals$`...`, args, check.attributes=FALSE))
  }
  stopifnot(all.equal(globals$sum, base::sum))
  stopifnot(all.equal(globals$x, x))

  message("\n*** globalsOf(dotdotdot='warn'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- globalsOf(expr, dotdotdot="warn")
  print(globals)
  if (name == "ok") {
    stopifnot(all.equal(names(globals), c("sum", "x")))
  } else {
    stopifnot(all.equal(names(globals), c("sum", "x", "...")))
    stopifnot(all.equal(globals$`...`, args, check.attributes=FALSE))
  }
  stopifnot(all.equal(globals$sum, base::sum))
  stopifnot(all.equal(globals$x, x))

  message("\n*** globalsOf(dotdotdot='error'):")
  cat(sprintf("Expression '%s':\n", name))
  print(expr)
  globals <- try(globalsOf(expr, dotdotdot="error"))
  if (name == "ok") {
    stopifnot(all.equal(names(globals), c("sum", "x")))
    stopifnot(all.equal(globals$sum, base::sum))
    stopifnot(all.equal(globals$x, x))
  } else {
    stopifnot(inherits(globals, "try-error"))
  }
} # for (name ...)

message("\n*** globalsOf(<exprs>, dotdotdot='return'):")
print(exprs)
globals <- globalsOf(exprs, dotdotdot="return")
print(globals)

} # aux()

aux(x=3:4, y=1, z=42L)
message("*** function(x, ...) globalsOf() ... DONE")


## Undo
options(opts)
