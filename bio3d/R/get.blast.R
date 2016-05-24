get.blast <- function(urlget, time.out = NULL, chain.single=TRUE) {

  if(substr(urlget, 1, 4) == "http" && grep("Blast.cgi", urlget)
      && grep("RID[[:space:]]*=", urlget)) {
     rid <- sub("^.*RID[[:space:]]*=[[:space:]]*", "", urlget)
     names(urlget)=rid 
  } else {
     stop("Illegal link for retrieving BLAST results")
  }

  cat(paste(" Searching ... please wait (updates every 5 seconds) RID =",rid,"\n "))

  ##- Retrieve results via RID code and check for job completion
  ##   (completion is based on retrieving HTML or CSV output)
  html <- 1
  t.count <- 0
  repeat {
    raw  <- try(read.csv(urlget,
                     header = FALSE, sep = ",", quote="\"", dec=".",
                     fill = TRUE, comment.char=""), silent=TRUE)
    if(class(raw)=="try-error") { stop("No hits found: thus no output generated") }
    html <- grep("DOCTYPE", raw[1,])
    
    if(!is.null(time.out) && (t.count > time.out) || (length(html) != 1))
      break;
    
    cat("."); Sys.sleep(5)
    t.count <- t.count + 5
  }
  
  if(length(html) == 1) {
     warning("\nTime out (", time.out, "s): Retrieve results with returned link\n", 
        urlget, "\n", sep="")
     return(urlget)
  }
  
  colnames(raw) <- c("queryid", "subjectids", "identity", "positives",
                     "alignmentlength", "mismatches", "gapopens",
                     "q.start", "q.end", "s.start", "s.end",
                     "evalue", "bitscore")
  
  ##- Expand 'raw' for each hit in 'subjectids' (i.e. split on ";")
  rawm <- as.matrix(raw)
  
  eachsubject <- strsplit(rawm[,"subjectids"],";")
  subjectids  <- unlist(eachsubject)
  n.subjects  <- sapply(eachsubject, length)
  
  rawm <- apply(rawm, 2, rep, times=n.subjects)
  rawm[,"subjectids"] <- subjectids
  

  ##- Parse ids
  all.ids <- strsplit(subjectids, "\\|")
  gi.id  <- sapply(all.ids, '[', 2)
  
  pdb.4char <- sapply(all.ids, '[', 4)
  pdb.chain <- sapply(all.ids, '[', 5)

  ## Catch long chain IDs as in hits from "P12612" (e.g "1WF4_GG" => "1WF4_g")
  if(chain.single) {
    chain.ind <- nchar(pdb.chain) > 1
    if(any(chain.ind)) {
      pdb.chain[ chain.ind ] <- tolower( substr(pdb.chain[ chain.ind ],1,1 ) )
    }
  }
  pdb.id <- paste(pdb.4char,"_",pdb.chain,sep="")


  ##- Map zero evalues to arbitrarily high value for -log(evalue)
  mlog.evalue <- -log(as.numeric(rawm[,"evalue"]))
  mlog.evalue[is.infinite(mlog.evalue)] <- -log(1e-308)

  
  cat(paste("\n Reporting",length(pdb.id),"hits\n"))
  
  output <- list(bitscore=  as.numeric(rawm[,"bitscore"]),
                evalue =  as.numeric(rawm[,"evalue"]),
                mlog.evalue = mlog.evalue,
                gi.id = gi.id,
                pdb.id = pdb.id,
                hit.tbl = rawm,
                raw = raw,
                url = urlget)

  class(output) <- "blast"
  return(output)
}
