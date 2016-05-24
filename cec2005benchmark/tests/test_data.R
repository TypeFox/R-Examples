library("cec2005benchmark")

.C("disablerand", PACKAGE = "cec2005benchmark")

extdata <- system.file("extdata", package = "cec2005benchmark")

for (i in 1:25) {
    test_data_file <- paste("test_data_func", i, ".txt", sep = "")
    test_data_path <- file.path(extdata, test_data_file)
    conn <- file(test_data_path, "r")
    X <- as.matrix(read.table(conn, nrows = 10, colClasses = "numeric"))
    f <- as.vector(as.matrix(read.table(conn, nrows = 10, colClasses = "numeric")))
    close(conn)
    fhat <- cec2005benchmark(i, X)
    stopifnot(isTRUE(all.equal(f, fhat, 1e-8)))
}
