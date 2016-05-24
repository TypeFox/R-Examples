
library(RcppCNPy)
library(rbenchmark)

## expensive: N <- 1e5
## cheaper:
n <- 1e4
k <- 50

M <- matrix(seq(1.0, n*k, by=1.0), n, k)

txtfile <- tempfile(fileext=".txt")
write.table(M, file=txtfile)

pyfile <- tempfile(fileext=".npy")
npySave(pyfile, M)

pygzfile <- tempfile(fileext=".npy.gz")
npySave(pygzfile, M)

print(do.call(rbind, (lapply(c(txtfile, pyfile, pygzfile),
                             function(f) file.info(f)["size"]))))

res <- benchmark(read.table(txtfile),
                 npyLoad(pyfile),
                 npyLoad(pygzfile),
                 order="relative",
                 columns=c("test", "replications", "elapsed", "relative"),
                 replications=10)
print(res)

unlink(txtfile)
unlink(pyfile)
unlink(pygzfile)



