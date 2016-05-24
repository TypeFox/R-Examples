Merge <- function(..., join="outer", fill=0, split=FALSE, verbose=TRUE)
{
   d <- list(...)
   cnames <- unique(unlist(lapply(d, colnames)))
   rnames <- unlist(lapply(d, rownames))
   dims <- lapply(d, dim)
   nd <- length(dims)
   nr <- sum(unlist(lapply(d, nrow)))
   METHODS <- c("inner", "outer", "leftouter")
   join <- pmatch(join, METHODS)
   if (is.na(join))
      stop(paste("Unknown join type: should either", paste("\"", METHODS, "\"", collapse=", ", sep="")))
   if (join==2) {
      x <- matrix(fill, nrow=nr, ncol=length(cnames)) 
      n <- 0
      for (i in 1:nd) {
         mt <- match(colnames(d[[i]]), cnames)
         nt <- (n+1):(n+dims[[i]][1])
         x[nt, mt] <- as.matrix(d[[i]])
         n <- n + dims[[i]][1]
      }
   } else if (join == 1) {
      allnames <- unlist(lapply(d, colnames))
      t <- table(allnames)
      cnames <- names(t)[t==nd]
      if (length(cnames)==0 & verbose)
         warning("One or more datasets have no varaibles in common.")
      cnames <- sort(cnames)
      x <- matrix(fill, nrow=nr, ncol=length(cnames))
      n <- 0
      for (i in 1:nd) {
        mt <- match(colnames(d[[i]]), cnames)
        mt2 <- na.omit(mt)
        if (length(mt2) == 0 & verbose)
           warning("One or more datasets have no variables in common.")
        nt <- (n+1):(n+dims[[i]][1])
        x[nt, mt2] <- as.matrix(d[[i]][, !is.na(mt)])
        n <- n + dims[[i]][1]
      }
   } else if (join == 3) {
     if (nd > 2) {
        d2 <- Merge(d[[-1]], fill=fill)
     } else {
        d2 <- d[[2]]
     }
     cnames <- sort(colnames(d[[1]]))
     ord <- order(colnames(d[[1]]))
     mt <- match(colnames(d2), cnames)
     mt2 <- na.omit(mt)
     if (length(mt2) == 0 & verbose)
        warning("One or more datasets have no variables in common.")
     x <- matrix(fill, nrow=nr, ncol=length(cnames))
     x[1:dims[[1]][1], ] <- as.matrix(d[[1]][, ord])
     x[(dims[[1]][1]+1):nr, mt2] <- as.matrix(d[[2]][, !is.na(mt)])
   }
   x <- as.data.frame(x)
   colnames(x) <- cnames
   rownames(x) <- make.unique(rnames)
   if (any(rownames(x) != rnames) & verbose)
      warning("Some row names were changed to avoid duplicates.")
   if (split) {
      l <- list(length(nd))
      n <- 0
      for (i in 1:nd) {
        nt <- (n+1):(n+dims[[i]][1])
        l[[i]] <- x[nt, ] 
        n <- n + dims[[i]][1]
      }
      x <- l
      args <- as.list(match.call())
      names(x) <- args[2:(1+nd)]
   }
   return(x)
}
