library("cec2005benchmark")

extdata <- system.file("extdata", package = "cec2005benchmark")

fbias_data_path <- file.path(extdata, "fbias_data.txt")
fbias <- as.vector(as.matrix(read.table(fbias_data_path, colClasses = "numeric")))
global_optima_path <- file.path(extdata, "global_optima.txt")
O <- as.matrix(read.table(global_optima_path, colClasses = "numeric"))

for (i in 1:25) {
    for (D in c(2, 10, 30, 50)) {
        o <- O[i, 1:D]
        if (i == 5) {
            o[1:ceiling(D/4)] <- -100 
            o[max(floor(0.75*D),1):D] <- 100
        }
        if (i == 8) {
            o[2*(1:floor(D/2)) - 1] <- -32
        }
        if (i == 20) {
            o[2*(1:floor(D/2))] <- 5
        }
        fhat <- cec2005benchmark(i, o)
        stopifnot(isTRUE(all.equal(fbias[i], fhat, 1e-8)))
    }
}
