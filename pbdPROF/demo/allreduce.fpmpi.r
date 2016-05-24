suppressMessages(library(pbdPROF, quietly = TRUE))

fn <- system.file("extdata/allreduce.fpmpi", package = "pbdPROF")
da <- read.prof(fn)

print(da)

