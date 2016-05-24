## from system.file("test-tools-1.R",      package = "Matrix"):
assertError <- function(expr) {
  d.expr <- deparse(substitute(expr))
  t.res <- tryCatch(expr, error = function(e) e)
  print(t.res)
  if(!inherits(t.res, "error"))
    stop(d.expr, "\n\t did not give an error", call. = FALSE)
  invisible(t.res)
}
assertWarning <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, warning = function(w)w)
    if(!inherits(t.res, "warning"))
	stop(d.expr, "\n\t did not give a warning", call. = FALSE)
    invisible(t.res)
}

library(HandTill2001)
data(ht01.twoclass)
new("bincap"
    , response = as.factor(ht01.twoclass$observed)
    , predicted = ht01.twoclass$predicted
    )


assertError(
  new("bincap"
      , response = ht01.twoclass$observed
      , predicted = ht01.twoclass$predicted
      )
  )
assertError(
  new("bincap"
      , response = as.factor(ht01.twoclass$observed)
      , predicted = as.data.frame(ht01.twoclass$predicted)
      )
  )
assertError(
  new("bincap"
      , true = 1
      )
  )
data(ht01.multipleclass)
new("multcap"
    , response = ht01.multipleclass$observed
    , predicted = as.matrix(ht01.multipleclass[, levels(ht01.multipleclass$observed)])
    )

assertError(
  new("multcap"
      , response = as.numeric(ht01.multipleclass$observed)
      , predicted = as.matrix(ht01.multipleclass[, levels(ht01.multipleclass$observed)])
      )
  )
assertError(
  new("multcap"
      , response = ht01.multipleclass$observed
      , predicted = ht01.multipleclass[, levels(ht01.multipleclass$observed)]
      )
  )
assertError(
  new("multcap"
      , true = 1
      )
  )
###########################################################

po <- ht01.multipleclass[, levels(ht01.multipleclass$observed)]
ro <- ht01.multipleclass$observed
######## get a  subset
i <- c(1,73,157,170,182,186)
p <- po[i,]
r <- ro[i]
multcap(
  response = r
  , predicted = as.matrix(p)
  )

## r extraneous levels
r <- ro[i]
p <- po[i,]
levels(r) <- c(levels(r), "foo", "bar")
assertWarning(
  multcap(
    response = r
    , predicted = as.matrix(p)
    )
)
## r levels unmatched by p
r <- ro[i]
p <- po[i,]
levels(r) <- c(levels(r), "foo", "bar")
assertWarning(
  multcap(
    response = r
    , predicted = as.matrix(p)
    )
)
## r values unmatched by p
r <- ro[i]
p <- po[i,1:5]
assertError(
  multcap(
    response = r
    , predicted = as.matrix(p)
    )
  )
## p columns unmatched by levels(r)
r <- ro[i]
p <- po[i,]
p$foo <- NA
assertError(
  multcap(
    response = r
    , predicted = as.matrix(p)
    )
  )
