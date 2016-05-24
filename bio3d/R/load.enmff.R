"ff.anm" <- function(r, cutoff=15, gamma=1, ...) {
  ifelse( r>cutoff, 0, gamma )
}

"ff.pfanm" <- function(r, cutoff=NULL, ...) {
  if(is.null(cutoff))
    return(r^(-2))
  else
    ifelse( r>cutoff, 0, r^(-2))
}

"ff.calpha" <- function(r, rmin=2.9, ...) {
  ## MMTK Units: kJ / mol / nm^2
  ##a <- 128; b <- 8.6 * 10^5; c <- 2.39 * 10^5;
  ## Bio3D Units: kJ / mol / A^2

  ## In case of unreasonable CA-CA distance
  if(!is.null(rmin))
    r[(r<rmin)] = rmin
  
  a <- 128 * 10^4; b <- 8.6 * 10^2; c <- 2.39 * 10^3;
  ifelse( r<4.0,
         b*(r) - c,
         a*(r)^(-6) )
}

"ff.calphax" <- function(r, atom.id, ssdat=NULL, verbose=FALSE, ...) {
  
  ## C-alpha FF (converted to kJ/mol/A**2)
  a <- 128 * 10^4
  b <- 8.6 * 10^2
  c <- 2.39 * 10^3
  
  ## Helix 1-4 interactions
  ah14.a <- -140.0318; ah14.b <- 922.1635;
  ah14.maxk <- 120;
  ah.func <- function(x) { ah14.a * x + ah14.b }
  
  ## Beta bridge interactions
  bb.a <- 49634261;
  bb.maxk <- 160.2831;
  bb.func <- function(x) { bb.a * x ^(-8) }
  
  ## SS-bonds
  ss.k <- 143.8896;
  
  ## Calculate default interactions
  ks <- ifelse( r<4.0, b * r - c, a * r^(-6))
  
  if(is.null(ssdat)) {
    return(ks)
  }
  else {
    ## type 1: helix 1-4 interactions
    tmp.inds <- which(ssdat$helix14[,1]==atom.id)
    inds.a   <- ssdat$helix14[tmp.inds, 2]
        
    ## type 2: beta bridges and strong hbonds
    tmp.inds <- which(ssdat$beta.bridges[,1]==atom.id)
    inds.b   <- ssdat$beta.bridges[tmp.inds, 2]
    
    ## type 3: ss-bonds
    tmp.inds <- which(ssdat$ss.bonds[,1]==atom.id)
    inds.c   <- ssdat$ss.bonds[tmp.inds, 2]
  }
  
  if ( length(inds.a) > 0 ) {
    dist.inds <- intersect(which(r<6.5), which(r>5.5))
    inds.a <- intersect(dist.inds, inds.a)
    
    if ( length(inds.a) > 0 ) {
      if(verbose) {
        cat("\n\nHelix 1-4 interactions:  ")
        cat(paste(atom.id, inds.a, sep=":"))
        cat("\n  Default values: \n   ")
        cat(paste(ks[inds.a], sep=" "))
      }
      
      ## Apply new force constants
      ks.tmp <- ah.func(r[inds.a])
      
      ## Avoid too large values
      ks.tmp[ which(ks.tmp>ah14.maxk) ] <- ah14.maxk
      
      ## Avoid values less than defaults
      old.vals <- ks[inds.a]
      reset.inds <- (ks.tmp < old.vals)
      ks.tmp[ reset.inds ] <- old.vals[ reset.inds ]
      ks[inds.a] <- ks.tmp
      
      if(verbose) {
        cat("\n  New interactions: \n   ")
        cat(ks[inds.a])
        cat("\n  Distances are: \n   ")
        cat(r[inds.a])
        cat("\n")
      }
    }
  }
  
  if ( length(inds.b) > 0 ) {
    dist.inds <- which(r<7)
    inds.b <- intersect(dist.inds, inds.b)
    
    if ( length(inds.b) > 0 ) {
      if(verbose) {
        cat("\n\nBeta bridge interactions:  ")
        cat(paste(atom.id, inds.b, sep=":"))
        cat("\n  Default values: \n   ")
        cat(ks[inds.b])
      }
      
      ## Apply new force constants
      ks.tmp <- bb.func(r[inds.b])
      
      ## Avoid too large values
      ks.tmp[ which(ks.tmp>bb.maxk) ] <- bb.maxk
      
      ## Avoid values less than defaults
      old.vals <- ks[inds.b]
      reset.inds <- (ks.tmp < old.vals)
      ks.tmp[ reset.inds ] <- old.vals[ reset.inds ]
      ks[inds.b] <- ks.tmp
      
      if(verbose) {
        cat("\n  New interactions: \n   ")
        cat(ks[inds.b])
        cat("\n  Distances are: \n   ")
        cat(r[inds.b])
        cat("\n")
      }
    }
  }
  
  if ( length(inds.c) > 0 ) {
    if(verbose) {
      cat("\n\nSS-bond:  ")
      cat(paste(atom.id, inds.c, sep=":"))
      cat("\n  Default values: \n   ")
      cat(ks[inds.c])
    }
    
    ks[inds.c] <- ss.k
    
    if(verbose) {
      cat("\n  New interactions: \n   ")
      cat(ks[inds.c])
      cat("\n  Distances are: \n   ")
      cat(r[inds.c])
      cat("\n")
    }
  }
  return(ks)
}


"ff.sdenm" <- function(r, atom.id, ssdat=NULL, ...) {
  ## sdENM by lazyload. contains an array with dimensions
  ## 20  x 20  x 27
  ## aa1 x aa2 x distance.category
  sdENM = bio3d::sdENM
 
  ## set sequence data to 1-letter aa code
  if(any(nchar(ssdat$seq)==3))
    sequ <- suppressWarnings( aa321(ssdat$seq) )  
  else
    sequ <- ssdat$seq

  ## Check for non-standard amino acids
  if(any(sequ=="X")) {
    cat("\n")
    unk <- paste(unique(ssdat$seq[which(sequ=="X")]), collapse=", ")
    stop(paste("Unknown aminoacid identifier for:", unk))
  }

  ## Initialize
  natoms <- length(r)
  aa.now <- sequ[atom.id]

  ## vector for spring constants
  ks <- rep(NA, natoms)
  
  ## Make distance categories
  map.dist <- c(0, seq(4, 16.5, by=0.5))
  
  dist.cat <- cut(r, breaks=map.dist, labels=FALSE)
  dist.cat[is.na(dist.cat)] <- 27
  dist.cat[atom.id]         <- NA

  ## Unique distance categories
  unq.cat <- unique(dist.cat)

  for ( i in 1:length(unq.cat) ) {
    tmp.cat <- unq.cat[i]
    
    if(!is.na(tmp.cat)) {
      tmp.inds <- which(dist.cat==tmp.cat)

      ## Since the lower.tri is NA we look up twice :P
      a <- sdENM[ aa.now, sequ, tmp.cat ][ tmp.inds ]
      b <- sdENM[ sequ, aa.now, tmp.cat ][ tmp.inds ]

      ks[ tmp.inds ][!is.na(a)] <- a[!is.na(a)]
      ks[ tmp.inds ][!is.na(b)] <- b[!is.na(b)]
    }
  }

  ## Set special restraints for covalent pairs
  inds.k12 <- c(atom.id -1, atom.id+1)
  inds.k12 <- inds.k12[ intersect(which(inds.k12 > 0), which(inds.k12 <= natoms)) ]
  ks[inds.k12] <- 43.52
  ks[atom.id]=0

  ## should in principle not get this far ...
  if(any(is.na(ks))) {
    stop(paste("Incompatible protein sequence:\n",
               " Paramters only exists for standard amino acid residues"))
  }
  
  ## sdENM FF is in arbitrary units
  ## The values given were arbitrarily normalized, so that
  ## the average kappa (over all amino acid pairs) is equal to 1, at d = 6 Ang.
  ## scale to kJ / mol / A^2 range:
  ks <- ks * 0.0083144621 * 300 * 10
  return(ks)
}

"ff.reach" <- function(r, atom.id, ssdat=NULL, ...) {
  natoms <- length(r)
  
  ## units in kJ/mol/A^2
  ## Table 1 - line DHFR
  #af <- 6770; as <- 2.08;
  #bf <- 0.951; bs <- 0.0589;
  #k12 <- 860; k13 <- 26.7; k14 <- 17;

  ## by correspondance with Kei (29 aug'13)
  ## line 38, page 1644, 2008 Biophysical J
  af <- 4810;  as <- 1.7;
  bf <- 0.872; bs <- 0.068;

  ## avgering over table 1
  k12 <- 866; k13 <- 28.7; k14 <- 24.16667;
  
  ## Calculate default interactions
  ks <- (af * exp(-bf*r)) + (as * exp(-bs*r))
  
  ## Differentiate between k12, k13, k14
  inds.k12 <- c(atom.id -1, atom.id+1)
  inds.k13 <- c(atom.id -2, atom.id+2)
  inds.k14 <- c(atom.id -3, atom.id+3)

  inds.k12 <- inds.k12[ intersect(which(inds.k12 > 0), which(inds.k12 <= natoms)) ]
  inds.k13 <- inds.k13[ intersect(which(inds.k13 > 0), which(inds.k13 <= natoms)) ]
  inds.k14 <- inds.k14[ intersect(which(inds.k14 > 0), which(inds.k14 <= natoms)) ]
  
  ks[inds.k12] <- k12;
  ks[inds.k13] <- k13;
  ks[inds.k14] <- k14;
  
  return(ks)
}


"load.enmff" <- function(ff='calpha') {
  
  ## Bahar "ANM"-ff
  if (ff=="anm")  {
    ff <- ff.anm
  }

  ## Yang Song and Jernigan (PNAS 2009)
  else if (ff=="pfanm") {
    ff <- ff.pfanm
  }
  
  ## Hinsen "C-alpha"-ff
  else if (ff=="calpha") {
    ff <- ff.calpha
  }
  
  ## Extended calpha-FF
  else if(ff=="calphax") {
    ff <- ff.calphax
  }

  ## sdENM
  else if(ff=="sdenm") {
    ff <- ff.sdenm
  }
  
  ## REACH
  else if(ff=="reach") {
    ff <- ff.reach
  }

  else {
    stop("force field not defined")
  }
  
  return(ff)
}
