# CHNOSZ/blast.R
# functions to analyze BLAST output files
# 20100320 jmd

## read a BLAST tabular output file, and filter by similarity,
## E-value and max hits per query
read.blast <- function(file, similarity=30, evalue=1e-5, max.hits=1, min.length=NA, quiet=FALSE) {
  # display some information about the file
  if("connection" %in% class(file)) fname <- summary(file)$description
  else fname <- basename(file)
  cat(paste("read.blast: reading", fname, "\n"))
  # read the blast tabular file
  blast <- read.csv(file,header=FALSE,sep="\t",stringsAsFactors=FALSE)
  if(!quiet) cat(paste("  read",nrow(blast),"lines with",length(unique(blast$V1)),"query sequences\n"))
  # take out rows that don't meet our desired similarity
  if(!is.na(similarity)) {
    is <- which(blast$V3 >= similarity)
    blast <- blast[is,]
    if(!quiet) cat(paste("  similarity filtering leaves",length(is),"lines and",length(unique(blast$V1)),"query sequences\n"))
  }
  # take out rows that don't meet our desired e-value
  if(!is.na(evalue)) {
    ie <- which(blast$V11 <= evalue)
    blast <- blast[ie,]
    if(!quiet) cat(paste("  evalue filtering leaves",length(ie),"lines and",length(unique(blast$V1)),"query sequences\n"))
  }
  # take out rows that don't have the minimum alignment length
  if(!is.na(min.length)) {
    ie <- which(blast$V4 >= min.length)
    blast <- blast[ie,]
    if(!quiet) cat(paste("  alignment length filtering leaves",length(ie),"lines and",length(unique(blast$V1)),"query sequences\n"))
  }
  # now take only max hits for each query sequence
  if(!is.na(max.hits)) {
    query.shift <- query <- blast$V1
    lq <- length(query)
    # for short (i.e., 1 query sequence) files, make sure that the hits get counted
    query.shift[max((lq-max.hits+1),1):lq] <- -1
    #query.shift <- query.shift[c((max.hits+1):lq,1:max.hits)]
    query.shift <- query.shift[c((lq-max.hits+1):lq,1:(lq-max.hits))]
    ib <- which(query!=query.shift)
    blast <- blast[ib,]
    if(!quiet) cat(paste("  max hits filtering leaves",length(ib),"lines and",length(unique(blast$V1)),"query sequences\n"))
  }
  # assign meaningful column names
  # http://bergmanlab.smith.man.ac.uk/?p=41, accessed 20111223
  colnames(blast) <- c("queryId", "subjectId", "percIdentity", "alnLength", "mismatchCount", 
    "gapOpenCount", "queryStart", "queryEnd", "subjectStart", "subjectEnd", "eVal", "bitScore")
  return(blast)
}

## process a blast table, identify the taxon
## for each hit
id.blast <- function(blast, gi.taxid, taxid.names, min.taxon=0, 
  min.query=0, min.phylum=0, take.first=TRUE) {
  # what are gi numbers of the hits
  # we use def2gi to extract just the gi numbers
  gi <- def2gi(blast[,2])
  # what taxid do they hit
  cat("id.blast: getting taxids ... ")
  imatch <- match(gi,gi.taxid[[1]])
  taxid <- gi.taxid[[2]][imatch]
  # what phyla are these
  cat("getting taxid.names ... ")
  iphy <- match(taxid,taxid.names$taxid)
  phylum <- taxid.names$phylum[iphy]
  species <- taxid.names$species[iphy]
  cat(paste(length(unique(phylum)),"unique phyla,",length(unique(taxid)),"unique taxa\n"))
  # now expand phyla into our blast table
  # we really don't want our stringsAsFactors (to use NAphylum, below)
  tax.out <- data.frame(taxid=taxid, phylum=as.character(phylum), species=species, 
    stringsAsFactors=FALSE)
  blast.out <- blast[,c(1,2,3,11)]
  colnames(blast.out) <- c("queryId","subjectId","percIdentity","eVal")
  blast <- cbind(blast.out,tax.out)
  # drop taxa that do not appear a certain number of times
  blast$taxid <- as.character(blast$taxid)
  nt <- table(blast$taxid)
  it <- which(nt/sum(nt) >= min.taxon)
  itt <- which(blast$taxid %in% names(nt)[it])
  blast <- blast[itt,]
  cat(paste("  min taxon abundance filtering leaves",length(unique(blast$queryId)),
    "query sequences,",length(unique(blast$phylum)),"phyla,",length(unique(blast$taxid)),"taxa\n"))
  # only take phylum assignments that make up at least a certain 
  # fraction ('amin') of hits to the query sequence
  uquery <- unique(blast$queryId)
  iquery <- match(uquery,blast$queryId)
  # function to select the (highest) represented phylum for each query
  iqfun <- function(i) {
    if((i-1)%%1000==0) cat(paste(i,""))
    myiq <- which(blast$query==uquery[i])
    # we don't have to do the calculation if there's only one hit
    if(length(myiq)==1) {
      iq <- iquery[i]
    } else {
      # take those hits, count each phyla, take the highest abundance,
      # check if it's above the minimum proportion of hits
      p <- as.character(blast$phylum[myiq])
      np <- table(p)
      pp <- np/sum(np)
      ip <- which.max(pp)
      if(pp[ip] < min.query) return(numeric())
      # use the first hit to that phylum
      myphy <- blast$phylum[myiq]
      iphy <- which(names(pp[ip])==myphy)
      iq <- myiq[iphy[1]]
    }
    return(iq)
  }
  # using lapply is tempting, but haven't got it working
  #iq <- as.numeric(lapply(1:length(uquery),iqfun))
  #iq <- iq[!is.na(iq)]
  if(min.query!=0) {
    cat("  filtering by phylum representation per query... ")
    iq <- numeric()
    for(i in 1:length(uquery)) iq <- c(iq,iqfun(i))
    cat("done!\n")
    blast <- blast[iq,]
    cat(paste("  min query representation filtering leaves",length(iq),"query sequences,",length(unique(blast$phylum)),
      "phyla,",length(unique(blast$taxid)),"taxa\n"))
  } else {
    # we'll just take the first hit for each query sequence
    if(take.first) blast <- blast[iquery,]
  }
  # now on to drop those phyla that are below a certain relative abundance
  # we count NA as a real phylum
  iNA <- which(is.na(blast$phylum))
  blast$phylum[iNA] <- "NAphylum"
  np <- table(blast$phylum)
  ip <- which(np/sum(np) >= min.phylum)
  ipp <- which(blast$phylum %in% names(np)[ip])
  blast <- blast[ipp,]
  cat(paste("  min phylum abundance filtering leaves",length(ipp),"query sequences,",length(unique(blast$phylum)),
    "phyla,",length(unique(blast$taxid)),"taxa\n"))
  return(blast)
} 

write.blast <- function(blast,outfile) {
  # using the output of read.blast
  # create a trimmed "blast format" output file with data only
  # in the following columns
  # 1 - query sequence ID
  # 2 - hit sequence ID (i.e., FASTA defline .. extract gi number)
  # 3 - similarity
  # 11 - evalue
  not <- rep(NA,length(blast[[1]]))
  not7 <- data.frame(not,not,not,not,not,not,not)
  c2 <- paste("gi",def2gi(blast[[2]]),sep="|")
  out <- data.frame(blast[[1]],c2,blast[[3]],not7,blast[[11]],not)
  write.table(out,outfile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE,na="")
  return(invisible(out))
}

def2gi <- function(def) {
  # extract gi numbers from FASTA deflines 20110131
  stuff <- strsplit(def,"\\|")
  gi <- sapply(1:length(stuff),function(x) {
    # the gi number should be in the 2nd position (after "gi")
    if(length(stuff[[x]])==1) return(stuff[[x]][1])
    else return(stuff[[x]][2])
  })
  return(gi)
}
