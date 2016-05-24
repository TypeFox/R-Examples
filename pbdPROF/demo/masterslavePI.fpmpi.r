suppressMessages(library(pbdPROF, quietly = TRUE))

fn <- system.file("extdata/masterslavePI.fpmpi", package = "pbdPROF")
da <- read.prof(fn)

print(da)

