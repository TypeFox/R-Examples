calcDistToRef <- function(megProj, spec, reference, align.exe, acc.tab, logfile, max.bp){
  
  conn <- dbconnect(megProj@db)
  
  obj <- as.list(dbReadDNA(conn, acc.tab, spec[1], max.bp,
                           enforce.binomial = FALSE))
  #         slog("\n--", length(obj), "sequences of", spec, 
  #             file = logfile)
  
  cat("\n--", length(obj), "sequences of", spec[1])
  
  ## select appropriate reference sequences
  ## --------------------------------------
  br <- spec[2]
  ## for small groups there might not be a benchamrk
  ## then choose the least sitant from availabel reference seqs
  if ( !br %in% names(reference) ){
    #slog("\n  - reference for", br, "not available ", file = logfile)
    obj <- c(reference, as.list(obj))
    obj <- mafft(obj, method = "auto", path = align.exe)
    d <- dist.dna(obj, model = "raw", pairwise.deletion = TRUE, 
                  as.matrix = TRUE)
    thisMean <- function(d, bn, sn){
      mean(d[bn, grep(sn, colnames(d))], na.rm = TRUE)
    }
    br <- sapply(names(reference), thisMean, d = d, sn = spec[1])
    br <- names(br)[which.min(br)]
    #slog("-- using", br, "instead", file = logfile)
  }
  obj <- c(reference[br], as.list(obj))
  obj <- mafft(obj, method = "auto", path = align.exe
               #                    , options = "--adjustdirectionaccurately")
  )
  ## delete positions that are 'gap' in reference
  obj <- obj[, obj[br, ] != as.raw(4)]
  
  ## calculate distance matrix
  d <- dist.dna(obj, model = "raw", pairwise.deletion = TRUE, 
                as.matrix = TRUE)
  d[is.na(d)] <- 1 ## assign distance = 1 to seqs that do not oberlap with reference seq
  string <- ""
  if ( nrow(d) == 2 & max(d, na.rm = TRUE) > .35 ){
    objRC <- revCompTest(obj, split = c(TRUE, FALSE), FALSE, align.exe)
    dd <- dist.dna(objRC, model = "raw", pairwise.deletion = TRUE, 
                   as.matrix = TRUE)
    if ( max(d) - max(dd) > .2 ) {
      d <- dd
      string <- " [RC]"
      dbWriteDNA(conn, acc.tab, del.gaps(objRC[grep(spec[1], rownames(objRC)), ]), 
                 enforce.binomial = FALSE)
    }
  }
  d <- d[, colnames(d) %in% br, drop = FALSE]
  d <- d[!rownames(d) %in% br, , drop = FALSE]
  
  ## write results to database
  ## -------------------------
  for ( j in rownames(d) ){
    
    gitax <- splitGiTaxon(j)
    gitax <- data.frame(gi = gitax["gi"], 
                        taxon = gitax["taxon"], 
                        distreference = min(d[j, ])
    )
    slog(paste("\n---- gi=" , gitax$gi, ": ", 
               round(gitax$distreference, 5), string, sep = ""), 
         file = logfile)
    gitax <- paste("UPDATE ", acc.tab, " SET distreference=", 
                   gitax$distreference, " WHERE gi='", 
                   gitax$gi, "' AND taxon='", 
                   gitax$taxon, "';", sep = "")
    dbSendQuery(conn, gitax)
  } # end of FOR-loop over j
  dbDisconnect(conn)
} 