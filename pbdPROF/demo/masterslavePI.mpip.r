suppressMessages(library(pbdPROF, quietly = TRUE))

fn <- system.file("extdata/masterslavePI.mpip", package = "pbdPROF")
da <- read.prof(fn)

print(da)

