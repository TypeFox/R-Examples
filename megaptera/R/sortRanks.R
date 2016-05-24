sortRanks <- function(x){
  
  ranks <- c("superkingdom", "kingdom",
             "phylum", "subphylum",
             "class", "subclass",
             "superorder", "order", "suborder",
             "superfamily", "family", "subfamily",
             "tribe", "subtribe",
             "genus", "subgenus",
             "species group",
             "species")
  
  nrank <- sapply(x, length)
  x <- x[nrank == max(nrank)]
  if ( length(x) == 1 ) return(x[[1]])
  
  ## these ranks have to be sorted
  r <- unique(unlist(x))
  r <- r[r != "no rank"]
  
  ## check for unknown ranks
  check <- r[!r %in% ranks]
  if ( length(check) > 0 ){
    stop("unknown rank: ", check)
  }
  
  ## sort ranks
  ranks <- ranks[ranks %in% r]
  r <- r[match(ranks, r)]

  ## insert internal "no rank" ranks
  ## -------------------------------
  y <- cbind(head(r, -1), tail(r, -1))
  id <- seq_along(r)
  nnr <- 0
  for ( k in 1:nrow(y) ){
    # index of ranks y[k, ] over elements in x (= rankSet)
    n <- sapply(x, match, x = y[k, ])
    n <- n[, !apply(n, 2, function(x) any(is.na(x)))]
    # remove duplicates
    n <- unique(t(n))
    # n := number of 'no rank'-ranks between y[k, ]
    if ( length(n)  == 0 ){
      n <- 0 # special case: no rankSet has both ranks
    } else {
      n <- max(n) - (min(n) + 1)
    }
    # nnr := cumulative sum of 'no rank'-ranks
    nnr <- nnr + n
    id[(k + 1):length(id)] <- id[(k + 1):length(id)] + n
  }
  rr <- rep("no rank", length(id) + nnr)
  rr[id] <- r
  
  ## leading "no rank" ranks
  ln <- max(sapply(x, function(x, y) which(x == y) - 1, y = rr[1] ))
  rr <- c(rep("no rank", ln), rr)
  rr
}



