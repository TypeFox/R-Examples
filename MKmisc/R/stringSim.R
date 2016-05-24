stringSim <- function (x, y, global = TRUE, match = 1, mismatch = -1, gap = -1,
                       minSim = 0){
    stopifnot(is.character(x), is.character(y))
    
    if (length(x) == 1 & nchar(x[1]) > 1)
        x1 <- strsplit(x, split = "")[[1]]
    else
        x1 <- x
    
    if (length(y) == 1 & nchar(y[1]) > 1)
        y1 <- strsplit(y, split = "")[[1]]
    else
        y1 <- y
    
    m <- length(x1)
    n <- length(y1)
    D <- matrix(NA, nrow = m+1, ncol = n+1)
    M <- matrix("", nrow = m+1, ncol = n+1)
    if(global){
      D[,1] <- seq_len(m+1)*gap+1
      D[1,] <- seq_len(n+1)*gap+1
      D[1,1] <- 0
    }else{
      D[,1] <- minSim
      D[1,] <- minSim
    }
    M[,1] <- "d"
    M[1,] <- "i"
    M[1,1] <- "start"
    if(global)
      text <- c("d", "m", "i")
    else
      text <- c("d", "m", "i", "stop")
    for(i in c(2:(m+1))){
      for(j in c(2:(n+1))){
        m1 <- D[i-1,j] + gap
        m2 <- D[i-1,j-1] + (x1[i-1] != y1[j-1])*mismatch + (x1[i-1] == y1[j-1])*match
        m3 <- D[i,j-1] + gap
        if(global){
          D[i,j] <- max(m1, m2, m3)
          wmax <- text[which(c(m1, m2, m3) == D[i,j])]
        }else{
          D[i,j] <- max(m1, m2, m3, minSim)
          wmax <- text[which(c(m1, m2, m3, minSim) == D[i,j])]
        }
        if("m" %in% wmax & x1[i-1] != y1[j-1])
          wmax[wmax == "m"] <- "mm"
        M[i,j] <- paste(wmax, collapse = "/")
      }
    }
    rownames(M) <- rownames(D) <- c("gap", x1)
    colnames(M) <- colnames(D) <- c("gap", y1)
    if(global)
      d <- D[m+1, n+1]
    else
      d <- max(D)
    
    attr(d, "Size") <- 2
    attr(d, "Diag") <- FALSE
    if(length(x) > 1) x <- paste0("", x, collapse = "")
    if(length(y) > 1) y <- paste0("", y, collapse = "")
    attr(d, "Labels") <- c(x,y)
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "ScoringMatrix") <- D
    attr(d, "TraceBackMatrix") <- M
    class(d) <- c("stringSim", "dist")

    return(d)
}
