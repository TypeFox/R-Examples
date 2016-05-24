require("slam")

x <- matrix(1 : 8, 2, 4)
dimnames(x) <- list(ROW = LETTERS[seq_len(nrow(x))],
                    COL = letters[seq_len(ncol(x))])
s <- as.simple_triplet_matrix(x)
dim(s) <- dim(x) <- c(4, 2)
stopifnot(identical(as.matrix(s), x))

d <- c(2, 3, 4)
x <- array(seq_len(prod(d)), d)
s <- as.simple_sparse_array(x)
dim(s) <- dim(x) <- c(d[length(d)], d[-length(d)])
stopifnot(identical(as.array(s), x))

dimnames(s) <- dimnames(x) <- NULL
stopifnot(identical(as.array(s), x))
