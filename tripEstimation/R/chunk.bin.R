`chunk.bin` <-
function (filename, pimgs, weights = NULL, chunk = 2000, proj = NULL)
{
    con <- file(filename, open = "rb")
    dm <- readBin(con, "integer", 3)
    cnt <- round(seq(0, 100, length = dm[3]/chunk))
    i <- 1
    cat("\n")
    repeat {
        n <- dm[1] * dm[2] * as.integer(chunk)
        A <- readBin(con, "double", n)
        m <- as.integer(length(A)/(dm[1] * dm[2]))
        if (m == 0)
            break
        A <- array(A, c(dm[-3], m))
	if (!is.null(proj)) A[,1:2,] <- apply(A[,1:2,], 3, function(x) project(x, proj))
	pimgs <- behav.bin(A, pimgs, weights = weights)
	cat(cnt[i], "...\n", sep = "")
	i <- i + 1
    }
    cat("\n")
    close(con)
    pimgs
}

