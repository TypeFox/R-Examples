print.mol2 <- function(x, ...) {
  cat(paste("... Name:", x$name, "\n"))
  cat(paste("...", nrow(x$atom), "atoms in molecule", "\n"))
  cat(paste("...", nrow(x$bond), "bonds in molecule", "\n"))
  cat(paste("...", nrow(x$xyz), "frame(s) stored", "\n"))
  cat("\n")

  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")
  cat("\n")
}


"read.mol2" <-
  function (file, maxlines = -1L)
  {
    if (missing(file)) {
      stop("read.mol2: please specify a MOL2 'file' for reading")
    }
    if (!is.numeric(maxlines)) {
      stop("read.mol2: 'maxlines' must be numeric")
    }
    toread <- file.exists(file)
    if (!toread) {
      stop("No input MOL2 file found: check filename")
    }
    
    atom.format <- matrix(c("eleno",   'numeric', 
                            "elena",   'character',
                            "x",       'numeric', 
                            "y",       'numeric', 
                            "z",       'numeric', 
                            "elety",   'character',
                            "resno",   'numeric', 
                            "resid",   'character',
                            "charge",  'numeric', 
                            "statbit", 'character'), ncol=2, byrow=TRUE,
                          dimnames = list(c(1:10), c("name","what")) )

    
    bond.format <-matrix( c("id",      'numeric', 
                            "origin",  'numeric', 
                            "target",  'numeric', 
                            "type",    'character',
                            "statbit", 'character'), ncol=2, byrow=TRUE,
                         dimnames = list(c(1:5), c("name","what")) )
    
    trim <- function(s) {
        s <- sub("^ +", "", s)
        s <- sub(" +$", "", s)
        s[(s == "")] <- NA
        s
      }

    
    split.line <- function(x, collapse=TRUE, ncol=NULL) {
      tmp <- unlist(strsplit(x, split=" "))
      inds <- which(tmp!="")
      if(!collapse)
        return(tmp[inds])
      else {
        tmp <- tmp[inds]
        if(length(tmp)<ncol)
          tmp <- c(tmp, NA)
        return(paste(tmp, collapse=";"))
      }
    }

    ## Read and parse mol2 file
    raw.lines <- readLines(file, n = maxlines)
    
    mol.start <- grep("@<TRIPOS>MOLECULE", raw.lines)
    atom.start <- grep("@<TRIPOS>ATOM", raw.lines)
    bond.start <- grep("@<TRIPOS>BOND", raw.lines)
    num.mol <- length(mol.start)

    if (!num.mol>0) {
        stop("read.mol2: mol2 file contains no molecules")
    }

    ## Fetch molecule names and info
    mol.names <- raw.lines[mol.start+1]
    mol.info <- trim( raw.lines[mol.start+2] )
    mol.info <- as.numeric(unlist(lapply(mol.info, split.line, collapse=FALSE)))

    ## mol.info should contain num_atoms, num_bonds, num_subs, num_feat, num_sets
    mol.info <- matrix(mol.info, nrow=num.mol, byrow=T)
 
    num.atoms <- as.numeric(mol.info[,1])
    num.bonds <- as.numeric(mol.info[,2])
    atom.end  <- atom.start + num.atoms
    bond.end  <- bond.start + num.bonds

    ## Build a list containing ATOM record indices
    se <- matrix(c(atom.start, atom.end), nrow=length(atom.start))
    atom.indices <- lapply(1:num.mol, function(d) seq(se[d,1]+1, se[d,2]))

    se <- matrix(c(bond.start, bond.end), nrow=length(bond.start))
    bond.indices <- lapply(1:num.mol, function(d) seq(se[d,1]+1, se[d,2]))

    ## Check if file consist of identical molecules
    same.mol <- TRUE
    mol.first <- NULL
    
    mols <- list()
    for ( i in 1:num.mol ) {
      raw.atom <- raw.lines[ atom.indices[[i]] ]
      raw.bond <- raw.lines[ bond.indices[[i]] ]
      
      ## Split by space
      txt <- unlist(lapply(raw.atom, split.line, ncol=10, collapse=TRUE))
      ncol <- length(unlist(strsplit(txt[1], ";")))
      atom <- read.table(text=txt, 
                         stringsAsFactors=FALSE, sep=";", quote='',
                         colClasses=atom.format[1:ncol,"what"],
                         col.names=atom.format[1:ncol,"name"],
                         comment.char="", na.strings=" ")


      ## Split by space
      txt <- unlist(lapply(raw.bond, split.line, ncol=5, collapse=TRUE))
      ncol <- length(unlist(strsplit(txt[1], ";")))
      bond <- read.table(text=txt,
                         stringsAsFactors=FALSE, sep=";", quote='',
                         colClasses=bond.format[1:ncol,"what"],
                         col.names=bond.format[1:ncol,"name"],
                         comment.char="", na.strings=" ")
      
      ## Same molecules as the previous ones?
      mol.str <- paste(atom$elena, collapse="")

      if ( i==1 ) {
        mol.first <- mol.str
      }
      else if (mol.str != mol.first) {
        same.mol <- FALSE
      }

      ## Store data
      xyz <- as.xyz(as.numeric(t(atom[, c("x", "y", "z")])))
      
      out <- list("atom" = atom, "bond" = bond, "xyz" = xyz, 
                  "info" = mol.info[i,], "name" = mol.names[i])
      class(out) <- "mol2"
      mols[[i]] <- out
    }

    ## If identical molecules
    if ( length(unique(num.atoms)) == 1 && same.mol == TRUE ) {
      atom <- mols[[1]]$atom
      bond <- mols[[1]]$bond
      xyz <- t(sapply(lapply(mols, function(x) x$xyz), rbind))
      xyz <- as.xyz(xyz)
      
      out <- list("atom" = atom, "bond" = bond, "xyz" = xyz,
                  "info" = mol.info[1,], "name" = mol.names[1])
      class(out) <- "mol2"
    }
    else {
      out <- mols
      ##class(out) <- "mol2s"
    }
    
    return(out)
  }
  
