# CHNOSZ/more.aa.R
# get amino acid compositions of proteins from model organisms
# (Eco.csv or Sce.csv)

more.aa <- function(protein=NULL, organism) {
  # return the composition of one or more proteins from
  # a "model organism", E. coli (Eco) or S. cerevisiae (Sce)
  # extracted from get.protein 20120519
  datapath <- paste("extdata/protein/", organism, ".csv.xz", sep="")
  datafile <- system.file(datapath, package="CHNOSZ")
  if(datafile=="") stop(paste("missing", datapath))
  mydata <- read.csv(datafile, as.is=TRUE)
  # if protein is not supplied, just give some information about the datafile
  if(is.null(protein)) {
    msgout("more.aa: ", datapath, " has data for ", nrow(mydata), " proteins\n")
    return(invisible())
  }
  # which columns to search for matches
  if(organism=="Sce") searchcols <- c("ORF", "SGDID", "GENE")
  else if(organism=="Eco") searchcols <- c("protein", "abbrv")
  # which columns have the amino acids, in the order of thermo$protein 
  iaa <- match(toupper(aminoacids(3)), toupper(colnames(mydata)))
  # iterate over a list
  waslist <- TRUE
  out <- list()
  if(!is.list(protein)) {
    waslist <- FALSE
    protein <- list(protein)
  }
  for(i in 1:length(protein)) {
    # find the matches
    imatch <- rep(NA, length(protein[[i]]))
    for(cname in searchcols) {
      icol <- match(cname, colnames(mydata))
      if(is.na(icol)) next
      iimatch <- match(protein[[i]], mydata[, icol])
      imatch[!is.na(iimatch)] <- iimatch[!is.na(iimatch)]
    }
    # report and remember the unsuccessful matches
    if(all(is.na(imatch))) stop("no proteins found!")
    inotmatch <- which(is.na(imatch)) 
    if(length(inotmatch) > 0) {
      if(length(inotmatch)==1) verb <- " was" else verb <- " were"
      msgout("more.aa: ", paste(protein[[i]][inotmatch], collapse=" "), verb, " not matched\n")
    }
    aa <- data.frame(mydata[imatch, iaa])
    # add the identifying columns
    if(organism=="Sce") ref <- mydata$SGDID[imatch]
    else ref <- rep(NA, length(protein[[i]]))
    if(organism=="Sce") abbrv <- mydata$GENE[imatch]
    else abbrv <- rep(NA, length(protein[[i]]))
    chains <- rep(1, length(protein[[i]]))
    chains[inotmatch] <- NA
    org <- rep(organism[[1]], length(protein[[i]]))
    precols <- data.frame(protein[[i]], organism=org, ref, abbrv, chains, stringsAsFactors=FALSE)
    colnames(precols)[1] <- "protein"
    colnames(aa) <- aminoacids(3)
    aa <- cbind(precols, aa)
    out <- c(out, list(aa))
  }
  # done!
  if(!waslist) return(out[[1]])
  else return(out)
}


