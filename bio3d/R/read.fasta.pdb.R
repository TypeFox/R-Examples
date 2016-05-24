"read.fasta.pdb" <-
function(aln, prefix="", pdbext="", fix.ali = FALSE, ncore=1, nseg.scale=1, ...) {

  ## Log the call
  cl <- match.call()
  
  # Parallelized by parallel package (Fri Apr 26 17:58:26 EDT 2013)
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  files  <- paste(prefix, aln$id, pdbext,sep="")

  ##cat(files,sep="\n")
  toread <- file.exists(files)

  ## check for online files
  toread[ substr(files,1,4)=="http" ] <- TRUE

  
  if(all(!toread))
    stop("No corresponding PDB files found")

  # Avoid multi-thread downloading
  if(any(substr(files,1,4) == "http")) {
     ncore = 1
     options(cores = ncore)
  }

  blank <- rep(NA, ncol(aln$ali))
 
  mylapply <- lapply
  if(ncore > 1) mylapply <- mclapply

#  for (i in 1:length(aln$id)) {
  retval <- mylapply(1:length(aln$id), function(i) {
    coords <- NULL; res.nu <- NULL
    res.bf <- NULL; res.ch <- NULL
    res.id <- NULL; res.ss <- NULL
    cat(paste("pdb/seq:",i,"  name:", aln$id[i]),"\n")

    if(!toread[i]) {
      warning(paste("No PDB file found for seq", aln$id[i],
              ": (with filename) ",files[i]), call.=FALSE)
      coords <- rbind(coords, rep(blank,3))
      res.nu <- rbind(res.nu, blank)
      res.bf <- rbind(res.bf, blank)
      res.ch <- rbind(res.ch, blank)
      res.id <- rbind(res.id, blank)
      res.ss <- rbind(res.ss, blank)
      
    } else {
      pdb <- read.pdb( files[i], verbose=FALSE, ... )
      ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
      pdbseq  <- aa321(pdb$atom$resid[ca.inds$atom])
      aliseq  <- toupper(aln$ali[i,])
      tomatch <- gsub("X","[A-Z]",aliseq[!is.gap(aliseq)])
      
      if(length(pdbseq)<1)
        stop(paste(basename(aln$id[i]), ": insufficent Calpha's in PDB"), call.=FALSE)
      
      ##-- Search for ali residues (1:15) in pdb
      start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])),collapse=""),
                           text = paste(pdbseq,collapse=""))[1]
            
      if (start.num == -1) {
        start.num <- 1
        ##stop(paste(basename(aln$id[i]), ": starting residues of sequence does not match starting residues in PDB"), call.=FALSE)
      }

      ##-- Numeric vec, 'nseq', for mapping aln to pdb
      nseq <- rep(NA,length(aliseq))
      ali.res.ind <- which(!is.gap(aliseq))
      if( length(ali.res.ind) > (length(pdbseq) - start.num + 1) ) {
        warning(paste(aln$id[i],
         ": sequence has more residues than PDB has Calpha's"), call.=FALSE)
        ali.res.ind <- ali.res.ind[1:(length(pdbseq)-start.num+1)] ## exclude extra
        tomatch <-  tomatch[1:(length(pdbseq)-start.num+1)]        ## terminal residues
      }
      nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))

      ##-- Check for miss-matches
      match <- aliseq != pdbseq[nseq] 
      if ( sum(match, na.rm=TRUE) >= 1 ) {
        mismatch.ind <- which(match)
        mismatch <- cbind(aliseq, pdbseq[nseq])[mismatch.ind,]
        n.miss <- length(mismatch.ind)

        if(sum(mismatch=="X") != n.miss) { ## ignore masked X res        
          details <- seqbind(aliseq, pdbseq[nseq])
          details$ali[is.na(details$ali)] <- "-"
          rownames(details$ali) <- c("aliseq","pdbseq")
          details$id <- c("aliseq","pdbseq")
          
          resmatch <- which(!apply(details$ali, 2, function(x) x[1]==x[2]))
          
          resid <- paste(pdb$atom$resid[ca.inds$atom][nseq][resmatch][1], "-",
                         pdb$atom$resno[ca.inds$atom][nseq][resmatch][1],
                         " (", pdb$atom$chain[ca.inds$atom][nseq][resmatch][1], ")", sep="")
          
          cat("\n ERROR   Alignment mismatch. See alignment below for further details\n")
          cat("         (row ", i, " of aln and sequence of '", aln$id[i], "').\n", sep="")
          cat("         First mismatch residue in PDB is:", resid, "\n")
          cat("         occurring at alignment position:", which(match)[1], "\n\n")
          .print.fasta.ali(details)
          
          msg <- paste(basename.pdb(aln$id[i]),
                       " alignment and PDB sequence miss-match\n",
                       "       beginning at position ",
                       which(match)[1], " (PDB RESNO ", resid, ")", sep="")
          stop(msg, call.=FALSE)
        }
      }
      
      ##-- Store nseq justified/aligned PDB data
      ca.ali <- pdb$atom[ca.inds$atom,][nseq,]
      coords <- rbind(coords, as.numeric( t(ca.ali[,c("x","y","z")]) ))
      res.nu <- rbind(res.nu, ca.ali[, "resno"])
      res.bf <- rbind(res.bf, as.numeric( ca.ali[,"b"] ))
      res.ch <- rbind(res.ch, ca.ali[, "chain"])
      res.id <- rbind(res.id, ca.ali[, "resid"])

      sse <- pdb2sse(pdb, verbose = FALSE)
      res.ss <- rbind(res.ss, sse[nseq])
    } ## end else for (non)missing PDB file
    return (list(coords=coords, res.nu=res.nu, res.bf=res.bf, res.ch=res.ch, res.id=res.id, res.ss=res.ss))
  } ) ## end mylapply

  retval <- do.call(rbind, retval)
  coords <- matrix(unlist(retval[, "coords"]), nrow=length(aln$id), byrow=TRUE)
  res.nu <- matrix(unlist(retval[, "res.nu"]), nrow=length(aln$id), byrow=TRUE)
  res.bf <- matrix(unlist(retval[, "res.bf"]), nrow=length(aln$id), byrow=TRUE)
  res.ch <- matrix(unlist(retval[, "res.ch"]), nrow=length(aln$id), byrow=TRUE)
  res.id <- matrix(unlist(retval[, "res.id"]), nrow=length(aln$id), byrow=TRUE)
 
  if( any(sapply(retval[, "res.ss"], is.null)) ) {
     res.ss <- NULL
  }
  else {
     res.ss <- matrix(unlist(retval[, "res.ss"]), nrow=length(aln$id), byrow=TRUE)
     rownames(res.ss) <- aln$id
  }

  rownames(aln$ali) <- aln$id
  rownames(coords) <- aln$id
  rownames(res.nu) <- aln$id
  rownames(res.bf) <- aln$id
  rownames(res.ch) <- aln$id
  rownames(res.id) <- aln$id
 
  if(fix.ali) {
     i1 <- which(is.na(res.nu))
     i2 <- which(is.gap(aln$ali))
     if(!identical(i1, i2)) {
        aln$ali[i1] <- aln$ali[i2[1]]
        warning("$ali component is modified to match $resno")
     }
  }

  out<-list(xyz=coords, resno=res.nu, b=res.bf,
            chain = res.ch, id=aln$id, ali=aln$ali, resid=res.id, sse=res.ss,
            call = cl)
  class(out)=c("pdbs","fasta")
  class(out$xyz) = "xyz"

  return(out)
}

