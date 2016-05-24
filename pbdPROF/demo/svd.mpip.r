suppressMessages(library(pbdPROF, quietly = TRUE))

fn <- system.file("extdata/svd.mpip", package = "pbdPROF")
da <- read.prof(fn)

print(da)

