"hmmer" <- function(seq, type='phmmer', db=NULL, verbose=TRUE, timeout=90) {
  cl <- match.call()
    
  oopsa <- requireNamespace("XML", quietly = TRUE)
  oopsb <- requireNamespace("RCurl", quietly = TRUE)
  if(!all(c(oopsa, oopsb)))
     stop("Please install the XML and RCurl package from CRAN")
 
  seqToStr <- function(seq) {
    if(inherits(seq, "fasta"))
      seq <- seq$ali
    if(is.matrix(seq)) {
      if(nrow(seq)>1)
        warning(paste("Alignment with multiple sequences detected. Using only the first sequence"))
      seq <- as.vector(seq[1,!is.gap(seq[1,])])
    }
    else
      seq <- as.vector(seq[!is.gap(seq)])
    return(paste(seq, collapse=""))
  }
  
  alnToStr <- function(seq) {
    if(!inherits(seq, "fasta"))
      stop("seq must be of type 'fasta'")
    tmpfile <- tempfile()
    write.fasta(seq, file=tmpfile)
    rawlines <- paste(readLines(tmpfile), collapse="\n")
    unlink(tmpfile)
    return(rawlines)
  }
  
  types.allowed <- c("phmmer", "hmmscan", "hmmsearch", "jackhmmer")
  if(! type%in%types.allowed )
    stop(paste("Input type should be either of:", paste(types.allowed, collapse=", ")))

  ## PHMMER (protein sequence vs protein sequence database)
  ## seq is a sequence
  if(type=="phmmer") {
    seq <- seqToStr(seq)
    if(is.null(db))
      db="pdb"
    db.allowed <- c("env_nr", "nr", "refseq", "pdb", "rp15", "rp35", "rp55",
                    "rp75", "swissprot", "unimes", "uniprotkb",
                    "uniprotrefprot", "pfamseq")
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop(paste("db must be either:", paste(db.allowed, collapse=", ")))
    
    seqdb <- db
    hmmdb <- NULL
    iter <- NULL
    rcurl <- TRUE
  }

  ## HMMSCAN (protein sequence vs profile-HMM database)
  ## seq is a sequence
  if(type=="hmmscan") {
    seq <- seqToStr(seq)
    if(is.null(db))
      db="pfam"
    db.allowed <- tolower(c("pfam", "gene3d", "superfamily", "tigrfam"))
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop(paste("db must be either:", paste(db.allowed, collapse=", ")))
    
    seqdb <- NULL
    hmmdb <- db
    iter <- NULL
    rcurl <- TRUE
  }

  ## HMMSEARCH (protein alignment/profile-HMM vs protein sequence database)
  ## seq is an alignment
  if(type=="hmmsearch") {
    if(!inherits(seq, "fasta"))
      stop("please provide 'seq' as a 'fasta' object")
    
    ##alnfile <- tempfile()
    ##seq <- write.fasta(seq, file=alnfile)
    seq <- alnToStr(seq)

    if(is.null(db))
      db="pdb"
    db.allowed <- tolower(c("pdb", "swissprot"))
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop(paste("db must be either:", paste(db.allowed, collapse=", ")))
    
    seqdb <- db
    hmmdb <- NULL
    iter <- NULL
    rcurl <- TRUE
  }

  ## JACKHMMER (iterative search vs protein sequence database)
  ## seq can be sequence, HMM, or multiple sequence alignment
  ## HMM not implemented here yet
  if(type=="jackhmmer") {
    if(!inherits(seq, "fasta"))
      stop("please provide 'seq' as a 'fasta' object")
    
    ##alnfile <- tempfile()
    ##seq <- write.fasta(seq, file=alnfile)
    seq <- alnToStr(seq)
    
    if(is.null(db))
      db="pdb"
    db.allowed <- tolower(c("pdb", "swissprot"))
    db.allowed <- c("env_nr", "nr", "refseq", "pdb", "rp15", "rp35", "rp55",
                    "rp75", "swissprot", "unimes", "uniprotkb",
                    "uniprotrefprot", "pfamseq")
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop(paste("db must be either:", paste(db.allowed, collapse=", ")))

    seqdb <- db
    hmmdb <- NULL
    iter <- NULL
    rcurl <- TRUE
  }
  
  
  ## Make the request to the HMMER website
  ##url <- paste('http://hmmer.janelia.org/search/', type, sep="")
  url <- paste("http://www.ebi.ac.uk/Tools/hmmer/search/", type, sep="")
  curl.opts <- list(httpheader = "Expect:",
                    httpheader = "Accept:text/xml",
                    verbose = verbose,
                    followlocation = TRUE
                    )
    
  hmm <- RCurl::postForm(url, hmmdb=hmmdb, seqdb=seqdb, seq=seq, 
                  style = "POST",
                  .opts = curl.opts,
                  .contentEncodeFun=RCurl::curlPercentEncode, .checkParams=TRUE )

  add.pdbs <- function(x, ...) {
    hit <- XML::xpathSApply(x, '@*')
    pdbs <- unique(XML::xpathSApply(x, 'pdbs', XML::xmlToList))
    new <- as.matrix(hit, ncol=1)
    
    if(length(pdbs) > 1) {
      for(i in 1:length(pdbs)) {
        hit["acc"]=pdbs[i]
        new=cbind(new, hit)
      }
      colnames(new)=NULL
    }
    return(new)
  }

  ##fetch.pdbs <- function(x) {
  ##  unique(XML::xpathSApply(x, 'pdbs', XML::xmlToList))
  ##}

  xml <- XML::xmlParse(hmm)
  data <- XML::xpathSApply(xml, '///hits', XML::xpathSApply, '@*')

  pdb.ids <- NULL
  if(db=="pdb") {
    tmp <- XML::xpathSApply(xml, '///hits', add.pdbs)
    data <- as.data.frame(tmp, stringsAsFactors=FALSE)
    colnames(data) <- NULL
  }

  data <- as.data.frame(t(data), stringsAsFactors=FALSE)
  data <- data[!duplicated(data$acc), ]
  ##rownames(data) <- data[, "acc"]
  
  ## convert to numeric
  fieldsToNumeric <- c("evalue", "pvalue", "score", "archScore", "ndom", "nincluded",
                       "niseqs", "nregions", "nreported", "bias", "dcl", "hindex")
  inds <- which(names(data) %in% fieldsToNumeric)
  
  for(i in 1:length(inds)) {
    tryCatch({
      data[[inds[i]]] = as.numeric(data[[inds[i]]])
    },
    warning = function(w) {
      #print(w)
      return(data[[inds[i]]])
    },
    error = function(e) {
      #print(e)
      return(data[[inds[i]]])
    }
    )
  }
  
  class(data) <- c("hmmer", type, "data.frame")
  out <- data
  return(out)
}
