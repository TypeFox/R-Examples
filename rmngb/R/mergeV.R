mergeV <- function(x, y,
                   by = intersect(names(x), names(y)), by.x = by, by.y = by,
                   all = FALSE, all.x = all, all.y = all,
                   verbose = TRUE, ...) {
  
  res <- merge(x = x, y = y,
               by.x = by.x, by.y = by.y,
               all.x = all.x, all.y = all.y, ...)
  if (verbose) {
    if (any(names(list(...)) == "incomparables"))
      warning("Not sure if it works when incomparable values are provided...")
    
    nX <- nrow(x)
    nY <- nrow(y)
    nR <- nrow(res)
    
    if (length(by.x) == 0L) {
      type <- "cross"
      tabCount <- data.frame(X = c(0, nX, 0, nX),
                             Y = c(0, nY, 0, nY),
                             R = c(0, nR, 0, nR),
                             row.names = c("X only", "X & Y", "Y only", "Total"))
    } else {
      type <- c("outer", "inner", "left", "right")[c(all.x && all.y, ! all.x && ! all.y,
                                                     all.x && ! all.y, ! all.x && all.y)]
      
      bx <- x[, by.x, drop = FALSE]
      by <- y[, by.y, drop = FALSE]
      
      names(bx) <- names(by) <- paste0("V", seq_len(ncol(bx)))
      
      bz <- do.call("paste", c(rbind(bx, by), sep = "\r"))
      bx <- bz[seq_len(nX)]
      by <- bz[nX + seq_len(nY)]
      
      if (type != "inner") {
        mR <- nrow(merge(x, y, by.x = by.x, by.y = by.y, all = FALSE, ...))
      } else {
        mR <- nR
      }
      tabCount <- data.frame(X = c(x1 <- sum(bx %out% by),
                                   sum(bx %in% by),
                                   0,
                                   nX),
                             Y = c(0,
                                   sum(by %in% bx),
                                   y1 <- sum(by %out% bx),
                                   nY),
                             R = c(if (type %in% c("outer", "left")) x1 else 0,
                                   mR,
                                   if (type %in% c("outer", "right")) y1 else 0,
                                   nR),
                             row.names = c("X only", "X & Y", "Y only", "Total"))
      
      tabX <- table(bx)
      tabY <- table(by)
      
      tabM <- merge(data.frame(id = names(tabX), X = as.vector(tabX)),
                    data.frame(id = names(tabY), Y = as.vector(tabY)),
                    all = FALSE)
      tabMatch <- with(tabM, table(X, Y))
    }
    tabCount[tabCount == 0] <- "."
    print(tabCount)
    cat("\nJoin type: ", type, "\n\n", sep = "")
    if (type != "cross")
      print(tabMatch)
  }
  res
}
