suppressMessages(library(pbdPROF, quietly = TRUE))

fn <- system.file("extdata/allreduce.mpip", package = "pbdPROF")
da <- read.prof(fn)

print(da)

