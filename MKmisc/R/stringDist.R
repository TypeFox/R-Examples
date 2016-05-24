stringDist <- function (x, y, method = "levenshtein", mismatch = 1, gap = 1){
    if (!is.na(pmatch(method, "levenshtein")))
        method <- "levenshtein"

    METHODS <- c("levenshtein", "hamming")
    method <- pmatch(method, METHODS)

    if (is.na(method))
        stop("invalid distance method")

    if (method == -1)
        stop("ambiguous distance method")

    stopifnot(is.character(x), is.character(y))
    
    if (length(x) == 1 & nchar(x[1]) > 1)
        x1 <- strsplit(x, split = "")[[1]]
    else
        x1 <- x
    
    if (length(y) == 1 & nchar(y[1]) > 1)
        y1 <- strsplit(y, split = "")[[1]]
    else
        y1 <- y
    
    if (method == 1){ ## Levenshtein
        m <- length(x1)
        n <- length(y1)
        D <- matrix(NA, nrow = m+1, ncol = n+1)
        M <- matrix("", nrow = m+1, ncol = n+1)
        D[,1] <- seq_len(m+1)*gap-1
        D[1,] <- seq_len(n+1)*gap-1
        D[1,1] <- 0
        M[,1] <- "d"
        M[1,] <- "i"
        M[1,1] <- "start"
        text <- c("d", "m", "i")
        for(i in c(2:(m+1))){
            for(j in c(2:(n+1))){
                m1 <- D[i-1,j] + gap
                m2 <- D[i-1,j-1] + (x1[i-1] != y1[j-1])*mismatch
                m3 <- D[i,j-1] + gap
                D[i,j] <- min(m1, m2, m3)
                wmin <- text[which(c(m1, m2, m3) == D[i,j])]
                if("m" %in% wmin & x1[i-1] != y1[j-1])
                  wmin[wmin == "m"] <- "mm"
                M[i,j] <- paste(wmin, collapse = "/")
            }
        }
        rownames(M) <- rownames(D) <- c("gap", x1)
        colnames(M) <- colnames(D) <- c("gap", y1)
        d <- D[m+1, n+1]
    }
    if(method == 2){ ## Hamming
        if(length(x1) != length(y1))
            stop("Hamming distance is only defined for equal length strings")
        d <- sum(x1 != y1)
        D <- NULL
        M <- NULL
    }
    attr(d, "Size") <- 2
    attr(d, "Diag") <- FALSE
    if(length(x) > 1) x <- paste0("", x, collapse = "")
    if(length(y) > 1) y <- paste0("", y, collapse = "")
    attr(d, "Labels") <- c(x,y)
    attr(d, "Upper") <- FALSE
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    attr(d, "ScoringMatrix") <- D
    attr(d, "TraceBackMatrix") <- M
    class(d) <- c("stringDist", "dist")

    return(d)
}
