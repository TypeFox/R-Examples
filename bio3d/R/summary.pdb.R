summary.pdb <- function(object, printseq=FALSE, ...) {

  ## Print a summary of basic PDB object features

  if( !is.pdb(object) ) {
    stop("Input should be a pdb object, as obtained from 'read.pdb()'")
  }

  ## Multi-model check and total atom count
  nmodel <- nrow(object$xyz)
  if( is.null(nmodel) ) {
    ntotal <- length(object$xyz)/3
    nmodel = 1
  } else {
    ntotal <- length(object$xyz[1,])/3
  }

  nxyz <- length(object$xyz)
  nres <- sum(object$calpha)
  chains <- unique(object$atom[,"chain"])

  all.inds <- atom.select(object, "all", verbose=FALSE)$atom
  prot.inds <- atom.select(object, "protein", verbose=FALSE)$atom
  nuc.inds <- atom.select(object, "nucleic", verbose=FALSE)$atom
  other.inds <- all.inds[! (all.inds %in% c(prot.inds, nuc.inds)) ]
  
  nprot <-length(prot.inds)
  nnuc <-length(nuc.inds)
  nresnuc <- length(unique(
    paste(object$atom$chain[nuc.inds], object$atom$insert[nuc.inds], object$atom$resno[nuc.inds], sep="-")))
                                                
  het <- object$atom[other.inds,]
  nhet.atom <- nrow(het)
  
  if(is.null(nhet.atom) | nhet.atom==0) {
    nhet.atom <- 0
    nhet.res <- 0
    hetres <- "none"
  } else { 
  	hetres.resno <- apply(het[,c("chain","resno","resid")], 1, paste, collapse=".")
  	nhet.res <- length(unique(hetres.resno))
	hetres.nres <- table(het[,c("resid")][!duplicated(hetres.resno)])
  	hetres <- paste( paste0( names(hetres.nres), " (",hetres.nres, ")"), collapse=", ")
  }

  if((nprot+nnuc+nhet.atom) != ntotal)
    warning("nPROTEIN + nNUCLEIC + nNON-PROTEIN + nNON-NUCLEIC != nTotal")

 
  cat("\n Call:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

  s <- paste0("\n   Total Models#: ", nmodel, 
  			  "\n     Total Atoms#: ", ntotal, ",  XYZs#: ", nxyz, 
  			  "  Chains#: ", length(chains),
  			  "  (values: ", paste(chains, collapse=" "),")",

             "\n\n     Protein Atoms#: ", nprot,
              "  (residues/Calpha atoms#: ", nres,")",
              
              "\n     Nucleic acid Atoms#: ", nnuc,
              "  (residues/phosphate atoms#: ", nresnuc,")",

             "\n\n     Non-protein/nucleic Atoms#: ", nhet.atom,
             "  (residues: ", nhet.res, ")",
             "\n     Non-protein/nucleic resid values: [ ", hetres," ]",
             "\n\n")
              
  cat(s)

  if(printseq) {
    ##protein
    if(nres>0) {
      prot.pdb <- trim.pdb(object, as.select(prot.inds))
      aa <- pdbseq(prot.pdb)
      if(!is.null(aa)) {
        if(nres > 225) {
          ## Trim long sequences before output
          aa <- c(aa[1:225], "...<cut>...", aa[(nres-3):nres])
        }
        aa <- paste("     ",  gsub(" ","", 
                                   strwrap( paste(aa,collapse=" "), 
                                           width=120, exdent=0) ), collapse="\n")
        cat("   Protein sequence:\n", aa, "\n\n", sep="")
      }
    }
    
    ## nucleic
    if(nresnuc>0) {
      na.pdb <- trim.pdb(object, as.select(nuc.inds))
      aa <- paste(object$atom$chain[nuc.inds], object$atom$insert[nuc.inds],
                  object$atom$resno[nuc.inds], object$atom$resid[nuc.inds],
                  sep="-")
      aa <- aa[!duplicated(aa)]
      aa <- unlist(lapply(strsplit(aa, "-"), function(x) x[4]))
      aa <- .aa321.na(aa)
      if(nresnuc > 225) {
        ## Trim long sequences before output
        aa <- c(aa[1:225], "...<cut>...", aa[(nresnuc-3):nresnuc])
      }
      aa <- paste("     ",  gsub(" ","", 
                                 strwrap( paste(aa,collapse=" "), 
                                         width=120, exdent=0) ), collapse="\n")
      cat("   Nucleic acid sequence:\n", aa, "\n\n", sep="")
    }
  }
    
  i <- paste( attributes(object)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")

  invisible( c(nmodel=nmodel, natom=ntotal, nxyz=nxyz, nchains=length(chains),  
               nprot=nprot, nprot.res=nres, nother=nhet.atom, nother.res=nhet.res) )

}
