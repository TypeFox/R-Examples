# stop()
stop("This is an error.\n")
a <- 3
if (a != 5) {stop("'a' should be equal to 5.\n")}
if (!(a == 5)) {stop("'a' should be equal to 5.\n")}
stopifnot(a == 5)

# warning()
warning("This is a warning.\n") 
if (a != 5) {warning("'a' should be equal to 5.\n")}
suppressWarnings(warning("This is a warning.\n"))

# message()
message("This is a diagnostic message.\n")
if (a != 5) {message("'a' should be equal to 5.\n")}
suppressMessages(message("This is a diagnostic message.\n"))

# return()
test.fun <- function(x) {
  y <- x + 10
  return(y)
}
test.fun(x = 23)

# options()
options(error = NULL)
options(error = dump.frames)
options(error = recover)

options(warn = -1)  # warnings ignored
options(warn = 0)   # default; warnings stored
options(warn = 1)   # warnings printed
options(warn = 2)   # Warnings turned into errors