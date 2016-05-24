"struct.aln" <-
  function(fixed, mobile, fixed.inds  = NULL, mobile.inds = NULL,
           write.pdbs = TRUE, outpath = "fitlsq",
           prefix = c("fixed", "mobile"),
           max.cycles = 10, cutoff = 0.5, ... ) {
    
    if(missing(fixed)) 
      stop("align: must supply 'pdb' object, i.e. from 'read.pdb'")
    if(missing(mobile)) 
      stop("align: must supply 'pdb' object, i.e. from 'read.pdb'")
    
    if(!is.pdb(fixed))
      stop("align: 'fixed' must be of type 'pdb'")
    if(!is.pdb(mobile))
      stop("align: 'mobile' must be of type 'pdb'")
    
    ## if indices are provided, make new PDB entities
    if ( !is.null(fixed.inds) ) {
      if(length(fixed.inds$atom)<2)
        stop("align: insufficent atom indices for fitting")
      
      #a <- NULL
      #a$atom <- fixed$atom[fixed.inds$atom, ]
      #a$xyz <- fixed$xyz[fixed.inds$xyz]
      #a$calpha <- as.logical(a$atom[,"elety"] == "CA")
      a <- trim.pdb(fixed, fixed.inds)
    }
    else {
      a <- fixed
      fixed.inds <- atom.select(fixed, 'all', verbose=FALSE)
    }
    
    if ( !is.null(mobile.inds) ) {
      if(length(mobile.inds$atom)<2)
        stop("align: insufficent atom indices for fitting")
      
      #b <- NULL
      #b$atom <- mobile$atom[mobile.inds$atom, ]
      #b$xyz <- mobile$xyz[mobile.inds$xyz]
      #b$calpha <- as.logical(b$atom[,"elety"] == "CA")
      b <- trim.pdb(mobile, mobile.inds)
    }
    else {
      b <- mobile
      mobile.inds <- atom.select(mobile, 'all', verbose=FALSE)
    }

    "xyz.dist" <- function(v) {
      a <- v[1:3]; b <- v[4:6]
      sqrt(sum((a-b)**2))
    }
    
    "resi.dev" <- function(xyz.a, xyz.b, cycle=1, cutoff = 0.5) {
      k <- matrix(xyz.a, ncol=3, byrow=T)
      l <- matrix(xyz.b, ncol=3, byrow=T)

      devs <- apply( cbind(k,l), 1, "xyz.dist")
      m <- median(devs)
      std <- sd(devs)
      
      cut <- m + (2*std)
      inds <- which( devs > cut )
      
      if ( (std < cutoff) || (length(inds)==0) ) {
        return( NULL )
      }
      else {
        cat( " Cycle ", i, ": ", length(inds), " atoms rejected", "\n", sep="")
        cat("  Mean: ", round(m,1),
            " Std: ", round(std,1),
            " Cut: ", round(cut,1), "\n", sep="" )
        return(inds)
      }
    }
    
    "remap.inds" <- function(pdb.init, inds.init, inds.trunc.atom) {
      ## Map back to indices for the entire PDB given
      inds.full <- NULL
      inds.full$atom <- inds.init$atom[inds.trunc.atom]
      inds.full$xyz <- atom2xyz(inds.full$atom)
      inds.full$logical <- atom2xyz(seq(1, nrow(pdb.init$atom))) %in% inds.full$xyz
      return(inds.full)
    }

    "parse.pdb" <- function(pdb, gaps, s, i) {
      pdbseq <- aa321(pdb$atom[pdb$calpha, "resid"])
      aliseq <- toupper(s$ali[i, ])
      tomatch <- gsub("X", "[A-Z]", aliseq[!is.gap(aliseq)])
      start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])), 
                             collapse = ""), text = paste(pdbseq, collapse = ""))[1]
      
      nseq <- rep(NA, length(aliseq))
      ali.res.ind <- which(!is.gap(aliseq))
      
      ali.res.ind <- ali.res.ind[1:length(pdbseq)]
      nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))
      
      pdb$atom <- cbind(pdb$atom, index=seq(1, nrow(pdb$atom)))
      ca.ali <- pdb$atom[pdb$calpha, ][nseq, ]
      at.inds <- ca.ali[, "index"]
      return(at.inds)
    }
    
    ## PDB list for sequence alignment
    pdb.list <- NULL
    pdb.list[[1]] <- a
    pdb.list[[2]] <- b
    
    ## Sequence alignment
    s <- lapply(pdb.list, pdbseq)
    s <- t(sapply(s, `[`, 1:max(sapply(s, length))))
    s[is.na(s)] <- "-"
    s <- seqaln(s, id = c("fixed", "mobile"), ...)
    gaps <- gap.inspect(s$ali)

    ## Parse truncated PDBs
    at.inds.a <- parse.pdb(a, gaps, s, 1)
    at.inds.b <- parse.pdb(b, gaps, s, 2)
    
    ## Fetch indices for fitting (truncated pdb)
    at.a <- as.numeric(at.inds.a[gaps$f.inds])
    at.b <- as.numeric(at.inds.b[gaps$f.inds])

    ## Indices for full pdb - done with the truncated ones
    a.inds.full <- remap.inds(fixed, fixed.inds, at.a)
    b.inds.full <- remap.inds(mobile, mobile.inds, at.b)

    ## Perform the initial fitting
    fit <- rot.lsq(mobile$xyz, fixed$xyz,
                   xfit=b.inds.full$logical, yfit=a.inds.full$logical)
    rmsd.init <- rmsd(as.vector(fixed$xyz), fit,
                      a.inds=a.inds.full$xyz, b.inds=b.inds.full$xyz)
    cat("\n")
    cat(" Initial RMSD (", length(gaps$f.inds), " atoms): ", rmsd.init,
        "\n", sep="")
    
    if ( write.pdbs ) {
      dir.create(outpath, FALSE)
      fname <- file.path(outpath, paste(prefix[2], "_", 0, ".pdb", sep=""))
      write.pdb(mobile, xyz=fit, file=fname)
    }

    ## Refinement process 
    rmsd.all <- c(rmsd.init)
    for ( i in seq(1,max.cycles) ) {
      if(i>max.cycles)
        break

      ## Find residues with largest structural deviation
      exc <- resi.dev(fixed$xyz[a.inds.full$xyz], fit[b.inds.full$xyz],
                      cycle = i, cutoff = cutoff)

      if ( is.null(exc) ) {
        break
      }
      else {
        ## Remove atoms for new round of fitting
        exc <- atom2xyz(exc)
        
        tmp <- seq(1,length( a.inds.full$logical ))
        exc.a <- tmp[which( a.inds.full$logical )][exc]
        a.inds.full$logical[exc.a] <- FALSE

        tmp <- seq(1,length( b.inds.full$logical ))
        exc.b <- tmp[which( b.inds.full$logical )][exc]
        b.inds.full$logical[exc.b] <- FALSE

        ## Build new xyz and atom indices
        a.inds.full$xyz <- which(a.inds.full$logical)
        b.inds.full$xyz <- which(b.inds.full$logical)

        a.inds.full$atom <- xyz2atom(a.inds.full$xyz)
        b.inds.full$atom <- xyz2atom(b.inds.full$xyz)
        
        ## Fit based on new indices
        fit <- rot.lsq(mobile$xyz, fixed$xyz,
                       xfit=b.inds.full$logical, yfit=a.inds.full$logical)

        if ( write.pdbs ) {
          fname <- file.path(outpath, paste(prefix[2], "_", i, ".pdb", sep=""))
          write.pdb(mobile, xyz=fit, file=fname)
        }
        
        ## Calculate RMSD 
        tmp.rmsd <- rmsd(as.vector(fixed$xyz), fit, a.inds=a.inds.full$xyz, b.inds.full$xyz)
        rmsd.all <- c(rmsd.all, tmp.rmsd)
        num.resi <- length(which(a.inds.full$logical))/3
        
        cat("  RMSD (", num.resi, " of ", length(gaps$f.inds), " atoms): ",
            tmp.rmsd, "\n", sep="")
      }
    }

    if ( write.pdbs ) {
      fname <- file.path(outpath, paste(prefix[1], ".pdb", sep=""))
      write.pdb(fixed, file=fname)
    }

    a.inds.full$logical <- NULL
    b.inds.full$logical <- NULL
    
    out <- list("a.inds"=a.inds.full, "b.inds"=b.inds.full,
                xyz=as.xyz(fit), rmsd=rmsd.all)
    
    return(out)     
  }
