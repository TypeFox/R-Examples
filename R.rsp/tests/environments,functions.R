library("R.rsp")

## Report template
report <- function() {
  cat("Report:\n")
  cat(sprintf("a=%g\n", a))
  cat("------\n")
} # report()


## 1. Generate report()
a <- 1
s <- rstring(report)
cat(s)
a <- 2
rcat(report)
a <- 3
dummy <- rfile(report, output="")


## 2. Generate report() in custom environment
env <- new.env()
env$a <- 1
s <- rstring(report, envir=env)
cat(s)

env$a <- 2
rcat(report, envir=env)

env$a <- 3
dummy <- rfile(report, output="", envir=env)



## 3. Generate report() within another function
foo <- function() {
  a <- 1
  s <- rstring(report)
  cat(s)
  a <- 2
  rcat(report)
  a <- 3
  dummy <- rfile(report, output="")
} # foo()

foo()
