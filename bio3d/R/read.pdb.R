"read.pdb" <-
function (file, maxlines=-1, multi=FALSE,
          rm.insert=FALSE, rm.alt=TRUE, ATOM.only = FALSE, verbose=TRUE) {

  if(missing(file)) {
    stop("read.pdb: please specify a PDB 'file' for reading")
  }
  if(!is.numeric(maxlines)) {
    stop("read.pdb: 'maxlines' must be numeric")
  }
  if(!is.logical(multi)) {
    stop("read.pdb: 'multi' must be logical TRUE/FALSE")
  }

  ##- Check if file exists locally or on-line
  toread <- file.exists(file)
  if(substr(file,1,4)=="http") { toread <- TRUE }

  ## Check for 4 letter code and possible on-line file
  if(!toread) {
    if(nchar(file)==4) {
      file <- get.pdb(file, URLonly=TRUE)
      cat("  Note: Accessing on-line PDB file\n")
    } else {
      stop("No input PDB file found: check filename")
    }
  }

  cl <- match.call()
  
  ## PDB FORMAT v3.3:    colpos, datatype,   name,       description
  atom.format <- matrix(c(6,    'character', "type",     # type(ATOM)
                          5,     'numeric',   "eleno",   # atom_no
                         -1,     NA,          NA,        # (blank)
                          4,     'character', "elety",   # atom_ty
                          1,     'character', "alt",     # alt_loc
                          4,     'character', "resid",   # res_na
                          1,     'character', "chain",   # chain_id
                          4,     'numeric',   "resno",   # res_no
                          1,     'character', "insert",  # ins_code
                         -3,     NA,           NA,       # (blank)
                          8,     'numeric',   "x",       # x
                          8,     'numeric',   "y",       # y
                          8,     'numeric',   "z",       # z
                          6,     'numeric',   "o",       # o
                          6,     'numeric',   "b",       # b
                         -6,     NA,           NA,       # (blank)
                          4,     'character', "segid",   # seg_id
                          2,     'character', "elesy",   # element_symbol
                          2,     'character', "charge"   # atom_charge (should be 'numeric']
                         ), ncol=3, byrow=TRUE,
                       dimnames = list(c(1:19), c("widths","what","name")) )

  trim <- function(s) {
    ##- Remove leading and trailing spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-""
    s
  }

  split.fields <- function(x) {
     ##- Split a character string for data.frame fwf reading
     ##  First splits a string 'x' according to 'first' and 'last'
     ##  then re-combines to new string with ";" as separator 
     x <- trim( substring(x, first, last) )
     paste(x,collapse=";")
  } 

  is.character0 <- function(x){length(x)==0 & is.character(x)}
 

  ##- Find first and last (substr) positions for each field
  widths <-  as.numeric(atom.format[,"widths"]) # fixed-width spec
  drop.ind <- (widths < 0) # cols to ignore (i.e. -ve)
  widths <- abs(widths)    # absolute vales for later
  st <- c(1, 1 + cumsum( widths ))
  first <- st[-length(st)][!drop.ind] # substr start
  last <- cumsum( widths )[!drop.ind] # substr end
  names(first) = na.omit(atom.format[,"name"])
  names(last) = names(first)


  ##- Read 'n' lines of PDB file
  raw.lines  <- readLines(file, n = maxlines)
  type <- substring(raw.lines, first["type"], last["type"])


  ##- Check number of END/ENDMDL records
  raw.end <- sort(c(which(type == "END"),
                    which(type == "ENDMDL")))

  ## Check if we want to store multiple model data
  if (length(raw.end) > 1) {
    cat("  PDB has multiple END/ENDMDL records \n")
    if (!multi) {
      cat("  multi=FALSE: taking first record only \n")
    } else {
      cat("  multi=TRUE: 'read.dcd/read.ncdf' will be quicker! \n")
      raw.lines.multi <- raw.lines
      type.multi <- type
    }
    raw.lines <- raw.lines[ (1:raw.end[1]) ]
    type <- type[ (1:raw.end[1]) ]
  }

  ##- Check for 'n' smaller than total lines in PDB file
  if ( length(raw.end) !=1 ) {
    if (length(raw.lines) == maxlines) {
      cat("  You may need to increase 'maxlines' \n")
      cat("   check you have all data in $atom \n")
    }
  }

  ##- Shortened records if ATOM.only = TRUE
  if(ATOM.only) {
     raw.lines <- raw.lines[type %in% c("HEADER", "ATOM  ", "HETATM")]
     type <- substring(raw.lines, first["type"], last["type"])
  }
 
  ##- Parse REMARK records for storing symmetry matrices to 
  ##  build biological unit by calling 'biounit()'
  remark <- .parse.pdb.remark350(raw.lines)

  ##- Split input lines by record type
  raw.header <- raw.lines[type == "HEADER"]
  raw.seqres <- raw.lines[type == "SEQRES"]
  raw.helix  <- raw.lines[type == "HELIX "]
  raw.sheet  <- raw.lines[type == "SHEET "]
  raw.atom   <- raw.lines[type %in% c("ATOM  ","HETATM")]

  if (verbose) {
    if (!is.character0(raw.header)) { cat(" ", raw.header, "\n") }
  }

  ## Edit from Baoqiang Cao <bqcao@ices.utexas.edu> Nov 29, 2009
  ## Old version:
  ##  seqres <- unlist(strsplit( trim(substring(raw.seqres,19,80))," +"))
  ## New version
  seqres <- unlist(strsplit( trim(substring(raw.seqres,19,80))," +"))
  if(!is.null(seqres)) {
    seqres.ch <- substring(raw.seqres, 12, 12)
    seqres.ln <- substring(raw.seqres, 13, 17)
    seqres.in <- ( !duplicated(seqres.ch) )
    names(seqres) <- rep(seqres.ch[seqres.in], times=seqres.ln[seqres.in])
  }
  ## End Edit from Baoqiang:

  ##- Secondary structure
  if(length(raw.helix) > 0) {
     helix  <- list(start = as.numeric(substring(raw.helix,22,25)),
                    end   = as.numeric(substring(raw.helix,34,37)),
                    chain = trim(substring(raw.helix,20,20)),
                    type  = trim(substring(raw.helix,39,40)))
   
     ##- insert code for initial and end residues of helices
     insert.i <- trim(substring(raw.helix,26,26))
     insert.e <- trim(substring(raw.helix,38,38))
     names(helix$start) <- insert.i
     names(helix$end) <- insert.e
  } else {
     helix <- NULL
  }
  if(length(raw.sheet) > 0) {
     sheet  <- list(start = as.numeric(substring(raw.sheet,23,26)),
                    end   = as.numeric(substring(raw.sheet,34,37)),
                    chain = trim(substring(raw.sheet,22,22)),
                    sense = trim(substring(raw.sheet,39,40)))
   
     ##- insert code for initial and end residues of sheets
     insert.i <- trim(substring(raw.sheet,27,27))
     insert.e <- trim(substring(raw.sheet,38,38))
     names(sheet$start) <- insert.i
     names(sheet$end) <- insert.e

     ##- remove repeated records for the same strand (e.g. in 1NH0)
     pa <- paste(sheet$start, insert.i, sheet$chain, sep='_')
     keep.inds <- which(!duplicated(pa))
     sheet <- lapply(sheet, '[', keep.inds)
  } else {
     sheet <- NULL
  }

 ## 2014-04-23
 ## Update to use single data.frame for ATOM and HETATM records
 ## file="2RGF"; multi=TRUE; 
 ## file="./4q21.pdb"; maxlines=-1; multi=FALSE; 
 ## rm.insert=FALSE; rm.alt=TRUE; het2atom=FALSE; verbose=TRUE

  atom <- read.table(text=sapply(raw.atom, split.fields), 
                    stringsAsFactors=FALSE, sep=";", quote='',
                    colClasses=atom.format[!drop.ind,"what"],
                    col.names=atom.format[!drop.ind,"name"],
                    comment.char="", na.strings="")

  ##-- End data.frame update


  ##- Coordinates only object
  ###xyz.models <- c(t(atom[,c("x","y","z")]))
  xyz.models <- matrix(as.numeric(c(t(atom[,c("x","y","z")]))), nrow=1)
  
  ##- Multi-model coordinate extraction 
  if (length(raw.end) > 1 && multi) {
      raw.atom  <- raw.lines.multi[ type.multi %in% c("ATOM  ","HETATM") ]

    if( (length(raw.atom)/length(raw.end)) ==nrow(atom) ){
      ## Only work with models with the same number of atoms)
      tmp.xyz=( rbind( substr(raw.atom, first["x"],last["x"]),
            substr(raw.atom, first["y"],last["y"]),
            substr(raw.atom, first["z"],last["z"]) ) )

      ## Extract coords to nrow/frame * ncol/xyz matrix
      xyz.models <- matrix( as.numeric(tmp.xyz), ncol=nrow(atom)*3, 
                            nrow=length(raw.end), byrow=TRUE)
      rownames(xyz.models) = NULL

    } else {
      warning(paste("Unequal number of atoms in multi-model records:", file))
    }
    rm(raw.lines.multi)
  }
  rm(raw.lines, raw.atom)


  ##- Possibly remove 'Alt records' (m[,"alt"] != NA)
  if (rm.alt) {
    if ( sum( !is.na(atom[,"alt"]) ) > 0 ) {
      first.alt <- sort( unique(na.omit(atom[,"alt"])) )[1]
      cat(paste("   PDB has ALT records, taking",first.alt,"only, rm.alt=TRUE\n"))
      alt.inds <- which( (atom[,"alt"] != first.alt) ) # take first alt only
      if(length(alt.inds)>0) {
        atom <- atom[-alt.inds,]
        xyz.models <- xyz.models[ ,-atom2xyz(alt.inds), drop=FALSE ] 
      }
    }
  }

  ##- Possibly remove 'Insert records'
  if (rm.insert) {
    if ( sum( !is.na(atom[,"insert"]) ) > 0 ) {
      cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
      insert.inds <- which(!is.na(atom[,"insert"])) # rm insert positions
      atom <- atom[-insert.inds,]
      xyz.models <- xyz.models[ ,-atom2xyz(insert.inds), drop=FALSE ]
    }
  }
  
  output<-list(atom=atom,
               #het=atom[atom$type=="HETATM",],
               helix=helix,
               sheet=sheet,
               seqres=seqres,
               xyz=as.xyz(xyz.models),
               calpha = NULL, remark = remark, call=cl)

  class(output) <- c("pdb", "sse")  
  ca.inds <-  atom.select.pdb(output, string="calpha", verbose=FALSE)
  output$calpha <- seq(1, nrow(atom)) %in% ca.inds$atom

  return(output)

}

##- parse REMARK records for building biological unit ('biounit()')
.parse.pdb.remark350 <- function(x) {

    raw.lines <- x

    # How many lines of REMARK 350?
    remark350 <- grep("^REMARK\\s+350", raw.lines)
    nlines <- length(remark350)

    # How many distinct biological unit?
    biolines <- grep("^REMARK\\s+350\\s+BIOMOLECULE", raw.lines)
    nbios <- length(biolines)

    if(nbios == 0) {
#       warning("REMARK 350 is incomplete.")
       return(NULL)
    }

    # End line number of each biological unit
    biolines2 <- c(biolines[-1], remark350[nlines])

    # How the biological unit was determined?
    method <- sapply(1:nbios, function(i) {
       author <- intersect(grep("^REMARK\\s+350\\s+AUTHOR DETERMINED BIOLOGICAL UNIT", raw.lines),
                            biolines[i]:biolines2[i])
       if(length(author) >= 1) return("AUTHOR")
       else return("SOFTWARE")
    } )
    # Get chain IDs to apply the transformation
    chain <- lapply(1:nbios, function(i) {
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines),
                            biolines[i]:biolines2[i])
       if(length(chlines) >= 1) {
          chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
          chs <- unlist(strsplit(chs, split=","))
       }
       else {
#          warning(paste("Can't determine chain IDs from REMARK 350 for biological unit",
#              i, sep=""))
          chs = NA
       }
       return(chs)
    } )
    if(any(is.na(chain))) return(NULL)

    mat <- lapply(1:nbios, function(i) {
       # Get transformation matrices
       mtlines <- intersect(grep("^REMARK\\s+350\\s+BIOMT", raw.lines),
                            biolines[i]:biolines2[i])
       # Get chain ID again: different trans matrices may be applied to different chains
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines),
                            biolines[i]:biolines2[i])
       chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
       chs <- strsplit(chs, split=",")

       if(length(mtlines) == 0 || length(mtlines) %% 3 != 0) {
#          warning("Incomplete transformation matrix")
          mat <- NA
       }
       else {
          mat <- lapply(seq(1, length(mtlines), 3), function(j) {
             mt <- matrix(NA, 3, 4)
             for(k in 1:3) {
                vals <- sub("^REMARK\\s+350\\s+BIOMT[123]\\s*", "", raw.lines[mtlines[j+k-1]])
                vals <- strsplit(vals, split="\\s+")[[1]]
                mt[k, ] <- as.numeric(vals[-1])
             }
             mt
          } )
          chs.pos <- findInterval(mtlines[seq(1, length(mtlines), 3)], chlines)
          names(mat) <- sapply(chs[chs.pos], paste, collapse=" ") ## apply each mat to specific chains
       }
       return(mat)
    } )
    if(any(is.na(mat))) return(NULL)

    out <- list(biomat = list(num=nbios, chain=chain, mat=mat, method=method))
    return(out)
}


