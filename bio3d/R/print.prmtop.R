print.prmtop <- function(x, printseq=TRUE, ...) {
  if(!is.null(x$SOLVENT_POINTER))
    sbox <- TRUE
  else
    sbox <- FALSE
  
  cn <- class(x)
  natom <-  x$POINTERS[1]
  ##nca <- length(which(x$ATOM_NAME=="CA" & x$AMBER_ATOM_TYPE=="CX"))
  
  if(sbox) {
    nres.total <- x$POINTERS[12]
    nres.solute <- x$SOLVENT_POINTER[1]

    nmol.total <- x$SOLVENT_POINTER[2]
    nmol.solute <- x$SOLVENT_POINTER[3]-1

    natom.per.mol <- x$ATOMS_PER_MOLECULE[1:nmol.solute]
    box.dim <- x$BOX_DIMENSIONS
  }
  else {
    nres.total <- x$POINTERS[12]
    nres.solute <- nres.total

    nmol.total <- length(which(x$ATOM_NAME=="OXT"))
    nmol.solute <- nmol.total
  }

  
  cat("\n Call:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  
  cat(" Class:\n  ", paste(cn, collapse=", "),
      "\n\n", sep = "")

  cat(" System information:", "\n")
  cat("      Total number of atoms:  ", natom, "\n", sep="")

  if(nres.solute!=nres.total) {
    cat("      Solute residues:  ", nres.solute, " (of ", nres.total, ")\n", sep="")
    cat("      Solute molecules:  ", nmol.solute, " (of ", nmol.total, ")\n", sep="")
  }
  else {
    cat("      Solute residues:  ", nres.solute, "\n", sep="")
    cat("      Solute molecules:  ", nmol.solute, "\n", sep="")
  }

  if(sbox)
    cat("      Box dimensions:  ", paste(round(box.dim,2), collapse=" x "), "\n", sep="")


  if(printseq) {
    aa <- aa321(x$RESIDUE_LABEL[1:nres.solute])
    if(nres.solute > 225) {
      ## Trim long sequences before output
      aa <- c(aa[1:225], "...<cut>...", aa[(nres.solute-3):nres.solute])
    }
    aa <- paste("     ",  gsub(" ","", 
            strwrap( paste(aa,collapse=" "), 
            width=120, exdent=0) ), collapse="\n")
    cat("\n")
    cat(" Sequence:\n", aa, "\n", sep="")


    ## other residues
    if(nres.total>nres.solute) {
      unq.res <- unique(x$RESIDUE_LABEL[(nres.solute+1):length(x$RESIDUE_LABEL)])
      unq.res <- paste("     ",  gsub(" ","", 
                       strwrap( paste(unq.res,collapse=" "), 
                        width=120, exdent=0) ), collapse="\n")
      cat("\n")
      cat(" Residues in solvent:\n", unq.res, "\n", sep="")
    }

  }

  

  cat("\n")

  #i <- paste( attributes(x)$names, collapse=", ")
  #cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")

  invisible( c(natom=natom, nres=nres.total, nres.solute=nres.solute,
               nmol=nmol.total, nmol.solute=nmol.solute) )
}
