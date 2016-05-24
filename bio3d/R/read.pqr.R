`read.pqr` <-
function (file, maxlines=-1, multi=FALSE,
          rm.insert=FALSE, rm.alt=TRUE, verbose=TRUE) {

  if(missing(file)) {
    stop("read.pqr: please specify a PQR 'file' for reading")
  }
  if(!is.numeric(maxlines)) {
    stop("read.pqr: 'maxlines' must be numeric")
  }
  if(!is.logical(multi)) {
    stop("read.pqr: 'multi' must be logical TRUE/FALSE")
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
                          8,     'numeric',   "o",       # o  ### 6 for pdb
                          8,     'numeric',   "b",       # b  ### 6 for pdb
                         -6,     NA,           NA,       # (blank)
                          4,     'character', "segid",   # seg_id
                          2,     'character', "elesy",   # element symbol
                          2,     'character', "charge"   # atom_charge (should be 'numeric']
                         ), ncol=3, byrow=TRUE,
                       dimnames = list(c(1:19), c("widths","what","name")) )

  trim <- function(s) {
    ##- Remove leading and trailing spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-NA
    s
  }

  split.fields <- function(x) {
     ##- Split a character string for data.frame fwf reading
     ##  First splits a string 'x' according to 'first' and 'last'
     ##  then re-combines to new string with "," as separator 
     x <- trim( substring(x, first, last) )
     paste(x,collapse=",")
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
    print("PDB has multiple END/ENDMDL records")
    if (!multi) {
      print("multi=FALSE: taking first record only")
    } else {
      print("multi=TRUE: 'read.dcd/read.ncdf' will be quicker!")
      raw.lines.multi <- raw.lines
      type.multi <- type
    }
    raw.lines <- raw.lines[ (1:raw.end[1]) ]
    type <- type[ (1:raw.end[1]) ]
  }

  ##- Check for 'n' smaller than total lines in PDB file
  if ( length(raw.end) !=1 ) {
    if (length(raw.lines) == maxlines) {
      print("You may need to increase 'maxlines'")
      print("check you have all data in $atom")
    }
  }

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
  helix  <- list(start = as.numeric(substring(raw.helix,22,25)),
                 end   = as.numeric(substring(raw.helix,34,37)),
                 chain = trim(substring(raw.helix,20,20)),
                 type  = trim(substring(raw.helix,39,40)))

  sheet  <- list(start = as.numeric(substring(raw.sheet,23,26)),
                 end   = as.numeric(substring(raw.sheet,34,37)),
                 chain = trim(substring(raw.sheet,22,22)),
                 sense = trim(substring(raw.sheet,39,40)))


 ## 2014-04-23
 ## Update to use single data.frame for ATOM and HETATM records
 ## file="2RGF"; multi=TRUE; 
 ## file="./4q21.pdb"; maxlines=-1; multi=FALSE; 
 ## rm.insert=FALSE; rm.alt=TRUE; het2atom=FALSE; verbose=TRUE

  atom <- read.table(text=sapply(raw.atom, split.fields),
                    stringsAsFactors=FALSE, sep=",", quote='',
                    colClasses=atom.format[!drop.ind,"what"],
                    col.names=atom.format[!drop.ind,"name"],
                    comment.char="")

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
  ##- Vector of Calpha positions 
  ##  check for calcium resid and restrict to ATOM records only
  calpha = (atom[,"elety"]=="CA") & (atom[,"resid"]!="CA") & (atom[,"type"]=="ATOM")

  output<-list(atom=atom,
               #het=atom[atom$type=="HETATM",],
               helix=helix,
               sheet=sheet,
               seqres=seqres,
               xyz=xyz.models,
               calpha = calpha, call=cl)

  class(output) <- c("pdb", "sse")
  class(output$xyz) <- c("numeric","xyz")
  return(output)

}
