"aa2mass" <-
  function(pdb, inds=NULL, mass.custom=NULL, addter=TRUE, mmtk=FALSE) {

    if (missing(pdb))
      stop("must supply 'pdb' object or vector of amino acid residue names")

    if(is.pdb(pdb)) {
      if(!is.null(inds)) {
        pdb <- trim.pdb(pdb, inds)
      }
      sequ <- pdb$atom[pdb$calpha,"resid"]
    }
    else {
      if(!is.null(inds))
        warning("'inds' has no effect when 'pdb' is vector")
      sequ <- pdb
      if(any(nchar(sequ)==1))
        sequ <- aa123(sequ)
      if(any(nchar(sequ)!=3))
        stop("must supply 'pdb' object or vector of amino acid residue names")
    }

    ## Define residues masses
    if(mmtk) {
      ## MMTK (for reproduction purposes!)
      w <- c( 71.079018, 157.196106, 114.104059, 114.080689, 103.143407,
             128.131048, 128.107678,  57.05203,  137.141527, 113.159985,
             113.159985, 129.18266,  131.197384, 147.177144,  97.117044,
              87.078323, 101.105312, 186.213917, 163.176449,  99.132996)
      
      aa <- c("ALA", "ARG", "ASN", "ASP", "CYS",
              "GLN", "GLU", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO",
              "SER", "THR", "TRP", "TYR", "VAL")

      mat <- data.frame(aa3=aa, aa1=aa321(aa), mass=w, formula=NA, name=NA)
      rownames(mat) <- aa
    }
    else  {
      ## Read data matrix
      mat <- bio3d::aa.table
    }
    
    ## Data frame with column names: aa3, aa1, mass, formula, name
    if (!is.null(mass.custom)) {
      if(class(mass.custom) != "list")
        stop("'mass.custom' must be of class 'list'")
      
      new.aas <- names(mass.custom)
      if(any(duplicated(new.aas))) {
        mass.custom[duplicated(new.aas)] <- NULL
        warning("duplicate residue name(s) in 'mass.custom'. using first occurrence(s) only.")
      }

      new.aas <- names(mass.custom)
      if(any(new.aas %in% mat$aa3)) {
        dups <- paste(unique(new.aas[new.aas %in% mat$aa3]), collapse=", ")
        warning(paste("residue name(s)", dups,
                      "exists in 'aa.table'. overwriting with provided value(s)."))
      }
      
      for(new.aa in new.aas) {
        if( new.aa %in% rownames(mat) ) {
          ## Replace residue mass
          mat[new.aa, "mass"] = mass.custom[[ new.aa ]]
        }
        else {
          ## Add new residue to data frame (aa.table)
          nr <- data.frame(list(aa3=new.aa, aa1="X",
                                mass=mass.custom[[ new.aa ]],
                                formula=NA, name=NA))
          rownames(nr) <- new.aa
          mat <- rbind(mat, nr)
        }
      }
    }

    ## Fetch mass from data frame
    wts <- mat[sequ, "mass"]
    
    ## Check for missing masses
    if(NA %in% wts) {
      inds <- which(wts %in% NA)
      unknown <- paste(unique(sequ[inds]), collapse=" ")
      stop(paste("Unknown amino acid identifier: ", unknown, sep=""))
    }
    
    if(addter) {
      wts[1] <- wts[1] + atom2mass("H")
      wts[length(wts)] <- wts[length(wts)] + atom2mass("O") + atom2mass("H")
    }
    return(wts)
  }
