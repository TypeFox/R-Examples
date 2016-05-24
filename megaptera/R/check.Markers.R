# LOCUS-WISE FREQUENCY OF SPECIES
# PACKAGE: megaptera
# CALLED BY: user
# AUTHOR: Christoph Heibl (at gmx.net)
# LAST CHANGE: 2014-11-01

check.Markers <- function(conn, subset, outgroup, colname = "blocks",
                          plot = TRUE){
    
    close.after <- FALSE
    if ( class(conn) == "megapteraProj" ){
      conn <- dbconnect(x@db)
      close.after <- TRUE
    }
    
    ## join taxonomy and locus tables
    ## ------------------------------
    x <- "SELECT * FROM taxonomy JOIN locus USING (spec)"
    x <- dbGetQuery(conn, x)
    
    if ( close.after ) dbDisconnect(conn)
    
    ## subset
    ## ------
    if ( !missing(subset) | !missing(outgroup) ){
      id.subset <- grep(paste(subset[-1], collapse = "|"), x[, subset[1]])
      if ( length(id.subset) == 0 ) stop("subset empty")
      id.outgroup <- grep(paste(outgroup[-1], collapse = "|"), x[, outgroup[1]])
      if ( length(id.outgroup) == 0 ) stop("outgroup empty")
      x <- x[ union(id.subset, id.outgroup), ]
    }
    colname <- paste("_", colname, sep = "")
    cols <- grep(colname, names(x))
    spec <- x$spec
    x <- x[, cols, drop = FALSE]
    
    ## some cosmetics on marker names
    colnames(x) <- gsub(colname, "", colnames(x))
    colnames(x) <- gsub("X_", "", colnames(x))
    colnames(x) <- gsub("([[:digit:]])(_)([[:digit:]])", "\\1,\\3", colnames(x))
    colnames(x) <- gsub("([[:alpha:]])(_)([[:alpha:]])", "\\1 \\3", colnames(x))
    bincov <- function(x){
      x[is.na(x)] <- 0
      x[grep("excluded", x)] <- 0
      x[grep("selected", x)] <- 1
      x[x > 1] <- 1
      as.numeric(x)
    }
    x <- apply(x, 2, bincov)
    rownames(x) <- spec
    
    ## number of species per marker
    sfreq <- sort(colSums(x), decreasing = TRUE)
      
    ## number of markers per species
    mfreq <- rowSums(x)
    private <- which(mfreq == 1)
    x <- x[private, , drop = FALSE]
    pfreq <- colSums(x)
    obj <- cbind(Ntotal = sfreq, Nprivate = pfreq[match(names(sfreq), names(pfreq))])
    
    pslist <- as.list(names(pfreq)[pfreq > 0])
    names(pslist) <- as.vector(pslist)
    for ( i in names(pslist) ){
      pslist[[i]] <- rownames(x)[x[, i] == 1]
    }
    
    ## produce barplot
    ## ---------------
    if ( plot ) {
      xx <- t(obj)
      xx[1, ] <- xx[1, ] - xx[2, ]
      opar <- par(no.readonly = TRUE)
      par(mar = c(4, 8, 4, 2))
      df.bar <- barplot(xx, horiz = TRUE, las = 1,
                        col = c("steelblue1", "orange"), border = NA,
                        main = "Number of species per marker",
                        legend.text = c("species represented by at least two markers", 
                                        "private (only this marker)"))
      id <- which(xx[2, ] > 0 )
      if ( length(id) > 0) {
        x <- (2 * xx[1, id] + xx[2, id]) / 2
        y <- df.bar[id]
        v <- xx[2, id]
        text(x, y, v)
      }
      par(opar)
    }
    
    ## return results
    ## --------------
    list(specPerMarker = obj, privateSpec = pslist)
}