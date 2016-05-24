mum <- function(x)
{
  if(!is.null(x)) {
    rn <- rownames(x)
    cn <- colnames(x)
    a <- duplicated(x[,1L], fromLast = FALSE)
    b <- duplicated(x[,1L], fromLast = TRUE)
    a[a != b] <- TRUE
    if(any(tot <- grepl(":total", rn))) {
      x <- x[!tot,]
      rn <- rn[!tot]
    }
    if(!is.matrix(x) && !is.null(x)) {
      x <- matrix(x, nrow = 1L)
      rownames(x) <- rn
      colnames(x) <- cn
    }
    if(nrow(x) > 1L) {
      a <- duplicated(x[,1L], fromLast = FALSE)
      b <- duplicated(x[,1L], fromLast = TRUE)
      a[a != b] <- TRUE
      if(any(a)) {
        drn <- rn[a]
        nc <- nchar(drn)
        rmrn <- drn[nc < max(nc, na.rm = TRUE)]
        x <- x[!(rn %in% rmrn),]
      }
      if(!is.null(x) && !is.matrix(x)) {
        x <- matrix(x, nrow = 1)
        rownames(x) <- rn[!(rn %in% rmrn)]
        colnames(x) <- cn
      }
    }
    rownames(x) <- rrmfs(rownames(x))
  }

  return(x)
}

