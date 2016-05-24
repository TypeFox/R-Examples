# PACKAGE: megaptera
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-11-01

stepG <- function(x, nob = FALSE){	
  
#   if ( !x$evaluate ) return(x)
  
  start <- Sys.time()
  quiet = FALSE
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## PARAMETERS
  ## ----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  spec.tab <- paste("spec", gsub("^_", "", gene), sep = "_")
  fract.miss <- x@params@fract.miss
  min.n.seq <- x@params@min.n.seq
  align.exe <- x@align.exe
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepG.log", sep = "-")
  if ( !quiet & file.exists(logfile) ) unlink(logfile)
  if ( !quiet )  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                      paste("\n", Sys.time(), sep = ""),
                      "\nSTEP G: alignment\n", file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
   
  ## clear <gene>_blocks attribute from previous runs of stepG/stepH
  ## ---------------------------------------------------------------
  SQL <- paste("UPDATE locus SET ", gene, "_blocks=NULL", sep = "")
  dbSendQuery(conn, SQL)
  
  ## read taxonomy relation from database
  ## ------------------------------------
  tax <- dbReadTaxonomy(conn, subset = spec.tab)
  
  ## alignment of genera -- either sequential or parallel
  ## ----------------------------------------------------
  gen <- unique(tax$gen)
  if ( length(gen) > 0 ) {
    cpus <- 12
    if ( length(gen) < cpus | !x@parallel ){
      seqs <- lapply(gen, alignGenus, x)
    } else {
      slog("\n", file = logfile)
      sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
      sfLibrary("megaptera", character.only = TRUE)
      megProj <- x
      sfExport("gen", "megProj")
      system.time(seqs <- sfLapply(gen, alignGenus, megProj = megProj))
      sfStop()
    }
  }
  names(seqs) <- gen
  
  ## if dataset contains only one genus, there is no need
  ## for alignment along guide tree
  ## ------------------------------
  if ( length(seqs) == 1 ){
    seqs <- seqs[[1]]
  } else {
    ## prepare guide tree for alignment
    ## --------------------------------
    gt <- tax2tree(tax, tip.rank = "gen")
    gt <- drop.tip(gt, pruned <- setdiff(gt$tip.label, names(seqs)))
    if ( length(pruned) > 0 ) 
      slog("\n\nPruning", length(pruned), "taxa from guide tree", file = logfile)
    seqs <- seqs[intersect(names(seqs), gt$tip.label)]
    
    ## loop over internal nodes of guide tree
    ## --------------------------------------
    slog("\n.. aligning sister clades higher than genera ..", 
                       file = logfile)
    i <- 1
    while ( Ntip(gt) > 2 ){
      #       for ( i in 1:7 ){
      #         load("~/r/dicaryota/data/BUGSEARCH.rda.RData")
      if ( !quiet ) slog("\n\nLEVEL", i, file = logfile)
      tp <- terminal.clades(gt)
      s1 <- seqs[gt$tip.label[-unlist(tp)]]
      names(seqs) <- match(names(seqs), gt$tip.label)
      ## align each set of sequences in list 'seqs'
      xx <- lapply(tp, function(tp, seqs) seqs[as.character(tp)], seqs = seqs)
      if ( !quiet ) slog(": aligning", length(xx), 
                         "pairs of sister clades", 
                         file = logfile)
      xx <- lapply(xx, mafft.merge, mafft.exe = align.exe)
      ## FILTER OUT WRONG SEQUENCES 
      xx <- lapply(xx, filter.alignment, megProj = x, logfile = logfile)
      
      if ( Nnode(gt) == 1 ) { ## basale Polytomie
        seqs <- xx
        i <- i + 1
        break
      }
      keep.tips <- sapply(tp, head, 1)
      drop.tips <- unlist(lapply(tp, tail, -1))
      gt$tip.label[keep.tips] <- 
        names(xx) <- paste("t-", i, "-", seq_along(xx), sep = "")
      gt <- drop.tip(gt, drop.tips)
      seqs <- c(xx, s1)
      
      if ( nob ){
        ## try to detect non-overlapping blocks
        md <- lapply(seqs, maxDistMPI)
        md <- which(md > .33)
        if ( length(md) > 0 ){
          if ( !quiet ) slog("\ntry to detect non-overlapping blocks", file = logfile)
          seqs[md] <- lapply(seqs[md], mafft, path = align.exe)
        }
      }
      
      i <- i + 1
    } # end of WHILE-loop
    if ( !quiet ) slog("\nLEVEL", i, file = logfile)
    #if ( Nnode(gt) > 1 ) # war das ein Fehler?
    # if ( Ntip(gt) > 1 ) 
    # REMARK: produziert einen Fehler wenn gt (aus der 
    # While-Schleife kommend) eine basale Polytomie hat: (A,
    # B, C, D);
    ## neuer Versuch:
    if ( length(seqs) > 2 )
      stop("uncomplete alignment in WHILE loop")
    if ( length(seqs) == 2 )
      seqs <- mafft.merge(seqs, align.exe)
    else seqs <- seqs[[1]]
  }
  
  ## prune ends of alignment from 'thin tails'
  ## -----------------------------------------
  seqs <- trimEnds(seqs, min(nrow(seqs), min.n.seq))
  seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  # cut positions from ends that have less than 3 unambigious bases
  #     n <- apply(seqs, 2, function(x) length(which(x %in% as.raw(c(136, 40, 72, 24)))))
  #     id <- which(n < 3)
  
  ## sort alignment taxonomically
  ## ----------------------------
  #seqs <- read.phy("_18s.phy")
  rownames(seqs) <- gsub("_R_", "", rownames(seqs))
  gt <- tax2tree(tax, tip.rank = "spec")
  gt <- ladderize(gt)
  gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(seqs)))
  seqs <- seqs[match(gt$tip.label, rownames(seqs)),]
  
  
  ## write alignments to PHY and NEX files
  ## -------------------------------------
  #     write.phy(seqs, paste(gene, "phy", sep = "."))
  #     rownames(seqs) <- gsub("-", "_", rownames(seqs))
  #     write.nex(seqs, paste(gene, "nex", sep = "."))
  write.dna.spectable(conn, spec.tab, seqs)
  ## This is just a hack to allow for neagative matching of regexs
  ## In the future stepG should insert the highest non-saturated taxonomic rank
  dbSendQuery(conn, paste("UPDATE", spec.tab, 
                          "SET block = 'included'", 
                          "WHERE block is NULL"))
  ## clear <gene>_blocks attribute from previous runs of stepG/stepH
  ## ---------------------------------------------------------------
  SQL <- paste("UPDATE locus ", 
               "SET ", gene, "_blocks='selected'",
               "WHERE ", sql.wrap(rownames(seqs), 
                                  "spec", "=", NULL),
               sep = "")
  lapply(SQL, dbSendQuery, conn = conn)
  dbDisconnect(conn)
  
  ## summary
  ## -------
  if ( !quiet ) {
    slog(
      paste("\n\nFinal alignment of", gene),
      paste("\nNumber of sequences:", nrow(seqs)),
      paste("\nNumber of base pairs:", ncol(seqs)),
      "\n\nSTEP G finished",
      file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  }
  
  invisible(x)
}
