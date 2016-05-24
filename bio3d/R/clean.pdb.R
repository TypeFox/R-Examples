#' Inspect And Clean Up A PDB Object
#'
#' Inspect alternative coordinates, chain breaks, bad residue
#' numbering, non-standard/unknow amino acids, etc. Return
#' a 'clean' pdb object with fixed residue numbering and optionally
#' relabeled chain IDs, corrected amino acid names, removed water,
#' ligand, or hydrogen atoms. All changes are recorded in a log in the 
#' returned object.
#'
#' @details call for its effects.
#' 
#' @param pdb an object of class \code{pdb} as obtained from
#'    function \code{\link{read.pdb}}. 
#' @param consecutive logical, if TRUE renumbering will result in
#'    consecutive residue numbers spanning all chains. Otherwise new residue
#'    numbers will begin at 1 for each chain.
#' @param force.renumber logical, if TRUE atom and residue records are renumbered 
#'    even if no 'insert' code is found in the \code{pdb} object.
#' @param fix.chain logical, if TRUE chains are relabeled based on chain breaks detected.
#' @param fix.aa logical, if TRUE non-standard amino acid names are converted into 
#'    equivalent standard names. 
#' @param rm.wat logical, if TRUE water atoms are removed.
#' @param rm.lig logical, if TRUE ligand atoms are removed.
#' @param rm.h logical, if TRUE hydrogen atoms are removed.
#' @param verbose logical, if TRUE details of the conversion process are printed.
#'
#' @return a 'pdb' object with an additional \code{$log} component storing 
#'    all the processing messages.
#'
#' @seealso \code{\link{read.pdb}}
#'
#' @author Xin-Qiu Yao & Barry Grant
#' 
#' @examples
#' \donttest{
#'    pdb <- read.pdb("1a7l")
#'    clean.pdb(pdb)
#' } 
"clean.pdb" <-
function(pdb, consecutive=TRUE, force.renumber = FALSE, fix.chain = FALSE, 
    fix.aa = FALSE, rm.wat = FALSE, rm.lig = FALSE, rm.h = FALSE, verbose=FALSE) {

  if(!is.pdb(pdb)) 
     stop("Input should be a 'pdb' object")

  cl <- match.call()

  ## processing message
  ## stored as an N-by-3 matrix with columns: 
  ##     FACT, OPERATION, IMPORTANT NOTE
  log <- NULL

  ## a flag to indicate if the pdb is clean
  clean <- TRUE

  ## Recognized amino acid names
  prot.aa <- bio3d::aa.table$aa3
  
  ## for residues and atoms renumbering
  first.eleno = 1
  first.resno = 1

  ## remove water
  if(rm.wat) {
     wat <- atom.select(pdb, "water", verbose = FALSE)
     if(length(wat$atom) > 0) {
        pdb$atom <- pdb$atom[-wat$atom, ,drop=FALSE]
        pdb$xyz <- pdb$xyz[, -wat$xyz, drop=FALSE]
        log <- .update.log(log, paste("Found", length(wat$atom), "water atoms"), "REMOVED")
     }
  }

  ## remove ligands 
  if(rm.lig) {
     lig <- atom.select(pdb, "ligand", verbose = FALSE)
     if(length(lig$atom) > 0) {
        pdb$atom <- pdb$atom[-lig$atom, ,drop=FALSE]
        pdb$xyz <- pdb$xyz[, -lig$xyz, drop=FALSE]
        log <- .update.log(log, paste("Found", length(lig$atom), "ligand atoms"), "REMOVED")
     }
  }

  ## remove hydrogens
  if(rm.h) {
     h.inds <- atom.select(pdb, "h", verbose = FALSE)
     if(length(h.inds$atom) > 0) {
        pdb$atom <- pdb$atom[-h.inds$atom, ,drop=FALSE]
        pdb$xyz <- pdb$xyz[, -h.inds$xyz, drop=FALSE]
        log <- .update.log(log, paste("Found", length(h.inds$atom), "hydrogens"), "REMOVED")
     }
  }
 
  ## check if 'alt' coords exist
  if(any(rm.p <- !is.na(pdb$atom$alt) & pdb$atom$alt != "A")) {
     pdb$atom <- pdb$atom[!rm.p, , drop=FALSE]
     pdb$xyz <- pdb$xyz[, -atom2xyz(which(rm.p)), drop=FALSE]
     log <- .update.log(log, paste("Found", sum(rm.p), "ALT records"), "REMOVED")
  }

  ## Some initial check: 
  ##   1. Are all amino acid and/or nucleic acid residues
  ##   distinguished by the combination chainID_resno_insert?
  .check.residue.ambiguity(pdb)

  ##   2. Check and clean up SSE annotation
  pdb <- .check.sse(pdb)
  log <- .update.log(log, pdb$log)

  ##   3. Fix pdb$calpha if it is mismatch pdb$atom
  ca.inds <- atom.select(pdb, "calpha", verbose = FALSE)
  calpha <- seq(1, nrow(pdb$atom)) %in% ca.inds$atom
  if(!identical(pdb$calpha, calpha)) {
     pdb$calpha <- calpha
     log <- .update.log(log, "pdb$calpha", "UPDATED")
  }
  
  ##   4. Fix object class if it is incorrect
  if(!inherits(pdb, "pdb") || !inherits(pdb, "sse") || 
       !inherits(pdb$xyz, "xyz")) {
     class(pdb) <- c("pdb", "sse")
     class(pdb$xyz) <- "xyz"
     log <- .update.log(log, "Object class", "UPDATED")
  }
  ###########
 
  ## following operations are on an independent object
  npdb <- pdb

  ## check chain breaks and missing chain ids
  has.fixed.chain <- FALSE
  capture.output( new.chain <- chain.pdb(npdb) )
  chain <- npdb$atom[, "chain"]  
  if(any(is.na(npdb$atom[, "chain"]))) {
     log <- .update.log(log, "Found empty chain IDs", 
               ifelse(fix.chain, "FIXED", "NO CHANGE"), 
               ifelse(fix.chain, "ALL CHAINS ARE RELABELED", ""))
     if(fix.chain) {
        npdb$atom[, "chain"] <- new.chain
        has.fixed.chain <- TRUE 
     } 
     else if(clean) 
        clean <- FALSE 
  }
  else {
     ## check if new chain id assignment is consistent to original one
     chn.brk <- bounds(chain[ca.inds$atom], dup.inds=TRUE, pre.sort=FALSE)
     new.chn.brk <- bounds(new.chain[ca.inds$atom], dup.inds=TRUE, pre.sort=FALSE)
     if(!isTRUE(all.equal(chn.brk, new.chn.brk))) {
        log <- .update.log(log, "Found inconsistent chain breaks",
                 ifelse(fix.chain, "FIXED", "NO CHANGE"),
                 ifelse(fix.chain, "ALL CHAINS ARE RELABELED", ""))
     
        log <- .update.log(log, "Original chain breaks:")
        if(nrow(chn.brk) == 1) {
           log <- .update.log(log, "  No chain break")
        } 
        else {
           pre.ca <- ca.inds$atom[chn.brk[-nrow(chn.brk), "end"]]
           pre.log <- capture.output( print(npdb$atom[pre.ca, c("resid", "resno", "chain")], row.names = FALSE) )
           log <- .update.log(log, pre.log)
        }
        log <- .update.log(log)
 
        log <- .update.log(log, "New chain breaks:")
        if(nrow(new.chn.brk) == 1) {
           log <- .update.log(log, "  No chain breaks")
        } 
        else { 
           new.ca <- ca.inds$atom[new.chn.brk[-nrow(new.chn.brk), "end"]]
           new.log <- capture.output( print(npdb$atom[new.ca, c("resid", "resno", "chain")], row.names = FALSE) )
           log <- .update.log(log, new.log)
        }

        if(fix.chain) {
           npdb$atom[, "chain"] <- new.chain
           has.fixed.chain <- TRUE 
        }
        else if(clean)
           clean <- FALSE
     } 
  }
   
  ## Renumber residues and atoms
  renumber <- FALSE
  if( any(!is.na(npdb$atom[, "insert"])) ) {
     renumber <- TRUE
     log <- .update.log(log, "Found INSERT records", "RENUMBERED")
     npdb$atom[, "insert"] <- as.character(NA)
  }
  else if(force.renumber) {
     renumber <- TRUE 
     log <- .update.log(log, "force.renumber = TRUE", "RENUMBERED")
  }
  else if(has.fixed.chain) {
     ## check again the ambiguity of residue labeling
     chk <- try(.check.residue.ambiguity(npdb))
     if(inherits(chk, "try-error")) {
        renumber <- TRUE
        log <- .update.log(log, "Found ambiguious residues after chain relabeling", "RENUMBERED")
     }
  }

  if(renumber) {

    ## Assign consecutive atom numbers 
    npdb$atom[,"eleno"] <- seq(first.eleno, length=nrow(npdb$atom))

    ## Determine what chain ID we have
    chain <- unique(npdb$atom[, "chain"])
    
    ##- Assign new (consecutive) residue numbers for each chain
    prev.chain.res = 0  ## Number of residues in previous chain
    for(i in 1:length(chain)) {
       inds <- which(npdb$atom[, "chain"] == chain[i])
       
       ## Combination of chain id, resno and insert code uniquely defines a residue (wwpdb.org)
       ## Here we use original pdb because we assume it should at least 
       ## distinguish different residues by above combination.
       ## We don't use the modified pdb (npdb) because all non-protein residues 
       ## are assigned a chain ID as "X" after calling chain.pdb(); 
       ## These residues could have the same resno (which are still in original 
       ## form) as they may be assigned different chain IDs in the original pdb. 
       res <- paste(pdb$atom[inds, "chain"], pdb$atom[inds, "resno"], 
                    pdb$atom[inds, "insert"], sep="_")

       n.chain.res <- length(unique(res))

       new.nums <- (first.resno+prev.chain.res):(first.resno+n.chain.res-1+prev.chain.res)
       npdb$atom[inds, "resno"] <- vec2resno(new.nums, res)

       if(consecutive) {
         ## Update prev.chain.res for next iteration 
         prev.chain.res = prev.chain.res + n.chain.res
       }
    }
  }

  ## update SSE
  if(has.fixed.chain || renumber) {
     ## Must use original pdb to unfold SSE
     sse <- pdb2sse(pdb, verbose = FALSE)  

     if(!is.null(sse)) {
        id <- sub(".*_.*_.*_([^_]*)$", "\\1", names(sse))
        names(sse) <- paste(npdb$atom[ca.inds$atom, "resno"], 
                            npdb$atom[ca.inds$atom, "chain"],
                            npdb$atom[ca.inds$atom, "insert"],
                            id, sep = "_")
        new.sse <- bounds.sse(sse)
     
        if(length(new.sse$helix$start) > 0) {
           npdb$helix$start <- new.sse$helix$start
           npdb$helix$end <- new.sse$helix$end
           npdb$helix$chain <- new.sse$helix$chain
        }
        if(length(new.sse$sheet$start) > 0) {
           npdb$sheet$start <- new.sse$sheet$start
           npdb$sheet$end <- new.sse$sheet$end
           npdb$sheet$chain <- new.sse$sheet$chain
        }
      
        if(!isTRUE(all.equal(npdb$helix, pdb$helix)) ||
           !isTRUE(all.equal(npdb$sheet, pdb$sheet)) )
           log <- .update.log(log, "SSE annotation", "UPDATED")
     }
  }
  
  ## update seqres 
  if(has.fixed.chain && !is.null(npdb$seqres)) {
     chs <- unique(npdb$atom[ca.inds$atom, "chain"])
     names(npdb$seqres) <- vec2resno(chs, names(npdb$seqres))
     if(!identical(npdb$seqres, pdb$seqres))
        log <- .update.log(log, "SEQRES", "UPDATED")
  }
 
  ## update amino acid name
  naa.atom <- which(npdb$atom[, "resid"] %in% prot.aa[-c(1:20)])
  naa.res <- intersect(ca.inds$atom, naa.atom)
  unk.atom <- which(!npdb$atom[, "resid"] %in% prot.aa)
  unk.res <- intersect(ca.inds$atom, unk.atom)
  if(length(naa.res) > 0) {
     log <- .update.log(log, paste("Found", length(naa.res), "non-standard amino acids"), 
              ifelse(fix.aa, "FIXED", "NO CHANGE"), 
              ifelse(fix.aa, "AMINO ACID NAMES ARE CHANGED", ""))
     tbl <- table(npdb$atom[naa.res, "resid"]) 
     tbl <- paste("  ", names(tbl), "(", tbl, ")", collapse = ",")
     log <- .update.log(log, tbl)
     if(fix.aa) {
         npdb$atom[naa.atom, "resid"] <- aa123(aa321(npdb$atom[naa.atom, "resid"])) 
         log <- .update.log(log, " Converted to")
         tbl <- table(npdb$atom[naa.res, "resid"])
         tbl <- paste("  ", aa123(aa321(names(tbl))), "(", tbl, ")", collapse = ",")
         log <- .update.log(log, tbl)
     }
     else 
        if(clean) clean <- FALSE
  }
  if(length(unk.res) > 0) {
     log <- .update.log(log, paste("Found", length(unk.res), "unknow amino acids"), 
               "NO CHANGE")
     tbl <- table(npdb$atom[unk.res, "resid"]) 
     tbl <- paste("  ", names(tbl), "(", tbl, ")", collapse = ",")
     log <- .update.log(log, tbl)
#     if(clean) clean <- FALSE
  } 

  ## update pdb$call
  npdb$call <- cl

  ## is the pdb clean?
#  if(clean)
#     log <- .update.log(log, "PDB is clean!")
  if(!clean) {
     msg <- "PDB is still not clean. Try fix.chain=TRUE and/or fix.aa=TRUE"
#     log <- .update.log(log, msg)
     warning(msg)
  }
     
  ## format log
#  log <- .format.log(log)
  if(verbose) print(log)

  npdb$log <- log
  return(npdb)
} 

.update.log <- function(log, fact="", op="", note="") {
  if(is.null(fact)) log
  else if(is.data.frame(fact)) .update.log(log, fact[,1], fact[,2], fact[,3])
  else rbind(log, data.frame(Data = fact, Action = op, Note = note))
}

.format.log <- function(log, format = c("print", "cat"), op.sign = "->", note.sign = "!!") {
  format <- match.arg(format)

  if(!is.null(log)) {
    log <- apply(log, 1, function(x) {
        if(nchar(x[2]) == 0)
           sprintf("%-40s", x[1])
        else if(nchar(x[3]) == 0)
           paste(sprintf("%-40s", x[1]), op.sign, sprintf("%-12s", x[2]))
        else
           paste(paste(sprintf("%-40s", x[1]), op.sign, sprintf("%-12s", x[2]), note.sign, x[3], note.sign) )
    } )
  }
  else {
    log <- "No problem found"
  }

  log <- switch(format, 
     print = log, 
     cat = paste(paste(log, collapse="\n"), "\n", sep="")
  )
  return(log)
}

.check.residue.ambiguity <- function(pdb) {

   ca.inds <- atom.select(pdb, "calpha", verbose = FALSE)
   c1p.inds <- atom.select(pdb, "nucleic", elety = "C1'", verbose = FALSE)
   inds <- combine.select(ca.inds, c1p.inds, operator = "+", verbose = FALSE)

   if(length(inds$atom) > 0) {
      strings <- paste(pdb$atom[inds$atom, "resno"],
                       pdb$atom[inds$atom, "chain"],
                       pdb$atom[inds$atom, "insert"], sep = "_")
      if(any(duplicated(strings)))
         stop(".check.residue.ambiguity(): Found ambiguous residue labeling")

   }
   invisible(NULL)
}

.check.sse <- function(pdb) {
   log <- NULL
   if(!is.null(pdb$helix) | !is.null(pdb$sheet)) {
      if(is.null(pdb$helix)) {
         pdb$helix <- list(start=NULL, end=NULL, chain=NULL, type=NULL) 
         log <- .update.log(log, "Helix is null but sheet is not", "UPDATED")
      }
      if(is.null(pdb$sheet)) {
         pdb$sheet <- list(start=NULL, end=NULL, chain=NULL, sense=NULL) 
         log <- .update.log(log, "Sheet is null but helix is not", "UPDATED")
      }
   }

   # if there is problem to generate sse vector
   ss <- try(pdb2sse(pdb, verbose = FALSE))
   if(inherits(ss, "try-error"))
      stop(".check.sse(): Unable to generate SSE sequence")

   pdb$log <- log

   return(pdb)
}
