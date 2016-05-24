barplot2D <-
function (z, x=0, y=0, width=1, height=1, colour, add=TRUE, col.frame=NULL, lwd.frame=1, threshold=1.1, ...) {
  area <- z
  stopifnot(is.vector(area), is.vector(colour),
            length(area) == length(colour),
            !all(is.na(area)))
  if (is.null(names(area))) {
    names(area) <- as.character(1:length(area))
  }
  area0 <- area
  if (any(is.na(area))) {
    warning("Discarding NAs")
    i <- which(!is.na(area))
    area <- area[i]
    colour <- colour[i]
  }
  stopifnot(all(area>=0), sum(area)>0)
  i <- order(-area)
  area <- area[i]
  colour <- colour[i]
  n <- length(area)
  res <- matrix(NA, nrow=n, ncol=8)
  colnames(res) <- as.vector(t(outer(LETTERS[1:4], 1:2, paste, sep="")))
  rownames(res) <- names(area)
  asp=1/(c(width,height)/max(width,height))
  x0=x-0.5*width
  y0=y-0.5*height
  x1=x+0.5*width
  y1=y+0.5*height
  A <- c(x0,y1)
  B <- c(x0,y0)
  C <- c(x1,y0)
  D <- c(x1,y1)
  if(add==F) {
    plot.new()
    plot.window(xlim=c(x0-width,x1+width), ylim=c(y0-height,y1+height)); axis(1); axis(2)
  }
  i <- 1
  while (i <= n) {
    lambda <- cumsum(area[i:n]) / sum(area[i:n])
    mu <- area[i]   / cumsum(area[i:n])
    nu <- area[i:n] / cumsum(area[i:n])
    penalty1 <- mu * sum(abs(A-B)*asp) / ( lambda * sum(abs(A-D)*asp) )
    penalty1 <- ifelse(penalty1 <= threshold, 0, penalty1 - threshold)
    penalty2 <- lambda * sum(abs(A-D)*asp) / ( nu * sum(abs(A-B)*asp) )
    penalty2 <- ifelse(penalty2 <= threshold, 0, penalty2 - threshold)
    j <- which.min(penalty1 + penalty2)[1] + i - 1
    lambda <- sum(area[i:j]) / sum(area[i:n])
    A1 <- A
    B1 <- B
    C1 <- (1-lambda) * B + lambda * C
    D1 <- (1-lambda) * A + lambda * D
    AA <- C1
    BB <- C
    CC <- D
    DD <- D1
    while (i <= j) {
      lambda <- area[i] / sum(area[i:j])
      B2 <- (1-lambda) * A1 + lambda * B1
      C2 <- (1-lambda) * D1 + lambda * C1
      polygon(rbind(A1, B2, C2, D1), col=colour[i], ...)
      res[i,] <- c(A1, B2, C2, D1)
      A1 <- B2
      D1 <- C2
      i <- i + 1
    }
    A <- AA
    B <- BB
    C <- CC
    D <- DD
  } # Main loop
  rect(x0,y0,x1,y1,border=col.frame,lwd=lwd.frame)
}

