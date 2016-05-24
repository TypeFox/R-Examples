# DETECTION and SEPARATION of unalignable SEQUENCE BLOCKS
# PACKAGE: megaptera
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-10-30

stepH <- function(x, clean = TRUE){
  
#   if ( !x$evaluate ) return(x)
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  spec.tab <- paste("spec", gsub("^_", "", gene), sep = "_")
  block.max.dist <- x@params@block.max.dist
  min.n.seq <- x@params@min.n.seq
  gblocks.exe <- x@mask.exe
  gb1 <- x@params@gb1
  gb2 <- x@params@gb2
  gb3 <- x@params@gb3
  gb4 <- x@params@gb4
  gb5 <- x@params@gb5
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepH.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP H: detection and separation of unalignable blocks\n", file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## read alignment
  ## --------------
  a <- dbReadDNA(conn, spec.tab, taxon = ".+", regex = TRUE,
                   ignore.excluded = TRUE)
  if ( !is.matrix(a) ){
    stop("stepG has not been called yet")
  }
  slog("\nalignment contains", nrow(a), "species", file = logfile)
  
  ## clear <gene>_blocks column from previous runs
  ## ---------------------------------------------
#   SQL <- paste("UPDATE locus SET ", gene, "_blocks=NULL", sep = "")
#   dbSendQuery(conn, SQL)
  
  ## FILTER OUT WRONG SEQUENCES 
  a <- list(filter.alignment(a, x, "G", logfile))
  
  ## custom distance function
  myDist <- function(d){
    if (!inherits(d, "DNAbin")) stop("object not of class \"DNAbin\"")
    n <- ncol(d)
    d <- dist.dna(d, model = "N", pairwise.deletion = TRUE, 
                  as.matrix = TRUE)
    d <- d/n
    d[d == 1] <- 0
    return(d)
  }
  
  ## detection of blocks
  ## -------------------
  d <- lapply(a, myDist)
  md <- sapply(d, max, na.rm = TRUE)
  md.id <- which(md >= block.max.dist)

  if ( length(md.id) == 0 ){
    slog("\n\nno blocks detected - nothing to do", file = logfile)
    SQL <- paste("UPDATE locus SET ", gene, "_blocks='selected' WHERE ", 
                 sql.wrap(rownames(a[[1]])), sep = "")
    dbSendQuery(conn, SQL)
    blk <- NULL
  } else {
    slog("\nat least 2 blocks detected", 
         "- starting separation of blocks:", 
         file = logfile)
    i <- 1
    slog("\n\tstep", i, "- maximum genetic distance in any block:", 
         max(md), file = logfile)
    while ( length(md.id) > 0 ){
#     for ( j in 1:8 ) {
      thisa <- a[[md.id[1]]]; thisd <- d[[md.id[1]]]
      ## species pair with largest genetic distance
      id <- which(thisd == max(thisd, na.rm = TRUE), arr.ind = TRUE)
      ## sometimes there are more than two species ...
      make2spec <- function(id){
        idd <- sort(table(rownames(id)))
        idd <- names(idd)[idd][1]
        id[c(which(id[, 2] == id[idd, "row"]), 
             which(rownames(id) == idd)), ]
      }
      while ( nrow(id) > 2) id <- make2spec(id)
      ## species that are genetically closer to first species of this pair
      id <- thisd[rownames(id)[1], ] < thisd[rownames(id)[2], ]
      id[is.na(id)] <- FALSE
      ## indices for both species
      id1 <- names(id)[id]; id2 <- names(id)[!id]
      a[[md.id[1]]] <- thisa[id1, ]
      a <- c(a, list(thisa[id2, ]))
      
      ## delete blocks containing less than <min.n.seq> species
      id <- sapply(a, nrow) >= min.n.seq
      if ( !all(id) ){
        SQL <- unlist(lapply(a[!id], rownames))
        SQL <- paste("UPDATE locus SET ", gene, "_blocks='excluded (stepH)' WHERE ", 
                     sql.wrap(SQL), sep = "")
        dbSendQuery(conn, SQL)
      }
      a <- a[id]
      if ( length(a) == 0 ) break
      d <- lapply(a, myDist)
      md <- sapply(d, max, na.rm = TRUE)
      md.id <- which(md >= block.max.dist)
      i <- i + 1
      slog("\n\tstep", i, "- maximum genetic distance in any block:", 
           max(md), file = logfile)
    } # end of WHILE
    a <- a[order(sapply(a, nrow), decreasing = TRUE)]
  } # end of ELSE
  
  if ( length(a) == 0 ){
    slog("\n\nall blocks smaller than <min.n.seq>", 
         "\nno alignment will be saved to file!",
         file = logfile)
  } else {
    ## cleaning of alignment
    ## ---------------------
    if ( clean ){
      slog("\n\nCleaning of sequence blocks\n", file = logfile)
      a <- lapply(a, gblocks, exec = gblocks.exe,
                  b1 = gb1, b2 = gb2, b3 = gb3, b4 = gb4, b5 = gb5) # with least conservative default
      a <- lapply(a, deleteEmptyCells, quiet = TRUE)
    }
    SQL <- lapply(a, rownames)
    SQL <- sapply(SQL, sql.wrap)
    SQL <- paste("UPDATE locus SET ", gene, "_blocks='selected (block-", 
                 1:length(SQL), ")' WHERE ", SQL, sep = "")
    lapply(SQL, dbSendQuery, conn = conn )
    
    blk <- sapply(a, nrow)
    
    ## concatenate alignment blocks
    ## ----------------------------
    if ( length(a) > 1 ) {
      # cbind.DNAbin should work now
      # a <- c.genes(single.list = a, match = FALSE)
      a <- do.call(cbind, c(a, fill.with.gap = TRUE))
    }
    else {
      a <- as.matrix(a[[1]])
      a <- deleteEmptyCells(a, quiet = TRUE)
      SQL <- paste("UPDATE", spec.tab, "SET", sql.wrap(0, term = "block"),
                   "WHERE", sql.wrap(rownames(a)))
      dbSendQuery(conn, SQL)
    }
    
    ## sort alignment taxonomically
    ## ----------------------------
    rownames(a) <- gsub("_R_", "", rownames(a))
    tax <- dbReadTaxonomy(conn, subset = a)
    gt <- tax2tree(tax, tip.rank = "spec")
    gt <- ladderize(gt)
    a <- a[match(gt$tip.label, rownames(a)), ]
    
    ## output of result
    ## ----------------
    slog(paste("\n\nFinal alignment of", gene),
         paste("\nNumber of sequences:", nrow(a)),
         paste("\nNumber of base pairs:", ncol(a)),
         file = logfile)
    if ( length(blk) > 0 ) slog("\nNumber of blocks:", length(blk), 
                                "(", paste(blk, collapse = ", "), ")", 
                                file = logfile)
    
    ## write to database
    ## -----------------
    write.dna.spectable(conn, spec.tab, a)
    
    ## write files
    ## -----------
    write.phy(a, paste(gene, "blocks.phy", sep = "."))
    rownames(a) <- gsub("-", "_", rownames(a))
    write.nex(a, paste(gene, "blocks.nex", sep = "."))
  }
  dbDisconnect(conn)
  slog("\n\nSTEP H finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  invisible(x)
}
