r2rmd <-
function(file, out=paste(file, "md", sep=""), header=TRUE, extra)
{
    trim <- 2
    if (missing(extra))
        extra <- ""
    x <- readLines(file)
    x2 <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
    n <- length(x)
    if (x2[n] != "") {
        x <- c(x, "")
        n <- n + 1
    }
    ch1 <- substr(x, 1, 1)
    ## dealing with less than trim leading hashes
    ch2 <- substr(x, 1, trim)
    ch3 <- unname(sapply(ch2, function(z) {
        val <- unique(sapply(1:nchar(z), function(i) substr(z, i, i)))
        if (length(val) > 1)
            z else val
        }))
    codecomment <- ch1 == "#" & ch3 != "#"
    ch1[codecomment] <- "c"
    ## find comments
    comment <- ch1 == "#"
    ## classify newlines
    if (ch1[1] %in% c("", "#"))
        comment[1] <- TRUE
    comment[n] <- TRUE
    for (i in 2:(n-1)) {
        if (comment[i-1] && ch1[i] %in% c("", "#"))
            comment[i] <- TRUE
    }
    for (i in (n-1):2) {
        if (comment[i+1] && ch1[i] %in% c("", "#"))
            comment[i] <- TRUE
    }
    ch4 <- substr(x, 1, trim + 1)
    trimchar <- c("## ", "##%", "###", "##")
    trim_comment <- comment & ch4 %in% trimchar
    #trim_comment[ch4 == paste0(rep("#", trim), collapse="")] <- ""
    if (header) {
        hv <- which(substr(x, trim+1, trim+3) == "---")[1:2]
        trim_comment[hv[1]:hv[2]] <- TRUE
    }
    ## chunks
    chunk <- integer(n)
    chunk[1] <- ifelse(comment[1], 0L, 1L)
    cmax <- chunk[1]
    for (i in 2:n) {
        if (comment[i-1] && !comment[i]) {
            cmax <- cmax + 1
            chunk[i] <- cmax
        }
        if (!comment[i-1] && !comment[i])
            chunk[i] <- cmax
    }
    nchunk <- max(chunk)
    ss <- matrix(0, nchunk, 2)
    for (i in 1:nchunk) {
        ss[i, 1] <- min(which(chunk == i))
        ss[i, 2] <- max(which(chunk == i))
    }
    x[ss[,2]] <- paste(x[ss[,2]], rep("\n```\n", nchunk), sep="")
    x[ss[,1]] <- paste(rep("```{r CHUNK_", nchunk),
        1:nchunk, rep(extra, nchunk),
        "}\n", x[ss[,1]], sep="")
    x[trim_comment] <- substr(x[trim_comment],
        trim + 1, nchar(x[trim_comment]))
    if (!is.null(out))
        writeLines(x, out)
    invisible(x)
}
