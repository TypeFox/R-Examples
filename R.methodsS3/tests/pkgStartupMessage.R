library("R.methodsS3")

message("TESTING: pkgStartupMessage()...")

msg <- "Hello world!"
pkgStartupMessage(msg)

for (quietly in c(NA, FALSE, TRUE)) {
  msg <- sprintf("Hello world! (quietly=%s)", quietly)
  pkgStartupMessage(msg, quietly=quietly)
}

message("TESTING: pkgStartupMessage()...DONE")
