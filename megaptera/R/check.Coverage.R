check.Coverage <- function(conn, pool.markers = TRUE){
  
  close.after <- FALSE
  if ( class(conn) == "sproj" ){
    conn <- dbConnect(PostgreSQL(), host = x$db[["host"]], 
                     port = x$db[["port"]], dbname = x$db[["dbname"]], 
                     user = x$db[["user"]], password = x$db[["password"]])
    close.after <- TRUE
  }
  
  x <- "SELECT * FROM taxonomy JOIN locus USING (spec)"
  x <- dbGetQuery(conn, x)
  if ( close.after ) dbDisconnect(conn)
  cols <- grep("_[gb|sel|blocks]", names(x))[-1]
  
  ## convert NA in _gb/_sel to 0
  ## and values > 0 to 1
  id <- cols[c(TRUE, TRUE, FALSE)]
  x[, id][is.na(x[, id])] <- 0
  x[, id][x[, id] > 0] <- 1
  
  ## coerce _blocks
  id <- cols[c(FALSE, FALSE, TRUE)]
  coerceBlocks <- function(b){
    id <- is.na(b) 
    if ( all(id) ){
      b[id] <- 1
    } else {
      id2 <- which(!id)
      id1 <- which(id)
      id3 <- grep("excluded", b[id2])
      if ( length(id3) > 0 ) id2 <- setdiff(id2, id3)
      b[union(id1, id3)] <- 0
      b[id2] <- 1
    }
    as.numeric(b)
  }
  x[, id] <- apply(x[, id], 2, coerceBlocks)
  
  if ( pool.markers ){
    tab <- data.frame(spec = x$spec, 
                      gb = rowSums(x[, cols[c(TRUE, FALSE, FALSE)]]),
                      sel = rowSums(x[, cols[c(FALSE, TRUE, FALSE)]]),
                      blocks = rowSums(x[, cols[c(FALSE, FALSE, TRUE)]])
    )
    tab$blocks[tab$blocks > tab$sel] <- tab$sel[tab$blocks > tab$sel]
    tab <- tab[order(tab$blocks, tab$sel, tab$gb, tab$spec), ]
    rownames(tab) <- NULL
  } else {
    id <- c(TRUE, FALSE, FALSE)
    tab <- data.frame(spec = x$spec, 
                      x[, cols[id]]
                      )
  }
  tab
}