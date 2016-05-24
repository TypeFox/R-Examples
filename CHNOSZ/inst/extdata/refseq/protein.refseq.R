# protein.refseq.R
# calculate the overall amino acid composition of proteins for each taxid in RefSeq
# 20100704 first version
# 20130922 deal with WP multispecies accessions (RefSeq61)

# uses system "join" command
# also, CHNOSZ for read.fasta, def2gi, etc.
require(CHNOSZ)

# where the microbial*.protein.faa.gz files are kept
faadir <- "protein"

# read WP multispecies table if it doesn't already exist
if(!exists("WP")) {
  cat("reading WP multispecies table\n")
  WP <- read.table(dir(pattern="multispecies_WP_accession_to_taxname.txt.gz"), sep="\t")
}

protein.refseq <- function(n=NULL) {
  # the list of sequence files (microbial44.protein.fa.gz etc)
  files <- dir(faadir)
  # loop over each one, getting the contents
  if(is.null(n)) n <- 1:length(files)
  out <- NULL
  for(i in n) {
    cat(paste("(", i, " of ", length(n), ") ", sep=""))
    # check if the name looks like a fasta file
    if(length(grep(".protein.faa.gz$", files[i]))==0) {
      cat(paste("skipping", files[i], "...\n"))
      next
    } else { 
      cat(paste("reading", files[i], "... "))
    }
    # read the file
    fa <- readLines(paste(faadir, files[i], sep="/"))
    # get the gi's from this file
    ihead <- grep("^>", fa)
    cat(paste("file has", length(ihead), "sequences\n"))
    # leave the gi's as character to sort lexigraphically (for use with system join command)
    gi <- def2gi(fa[ihead])
    # now find the taxids that belong to these gi's
    cat("finding matching taxids in catalog ...")
    write.table(sort(gi), "gi.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
    system("join -j 1 gi.txt gi.taxid.txt > gi.taxid.match")
    gi.taxid <- read.table("gi.taxid.match", header=FALSE)
    # produce an expanded table, replacing multispecies WP_ records with all corresponding taxids
    is.multi <- gi.taxid$V1 %in% WP$V2
    is.in.gi <- WP$V2 %in% gi.taxid$V1
    # take out generic records
    gi.taxid <- gi.taxid[!is.multi, ]
    # add specific records
    gi.taxid.WP <- WP[is.in.gi, ]
    colnames(gi.taxid.WP) <- c("V0", "V1", "V2", "V3")
    gi.taxid <- rbind(gi.taxid, gi.taxid.WP[, c(2, 3)])
    # identify unique taxids
    utax <- unique(gi.taxid$V2)
    cat(paste(" found", length(utax), "unique taxa\n"))
    # now loop over each taxon
    for(j in 1:length(utax)) {
      # report on the taxon name and number of sequences
      jtax <- which(gi.taxid$V2==utax[j])
      nseq <- length(jtax)
      cat(paste("taxon", utax[j], "has", nseq, "sequences "))
      # match the gi numbers to their position in the fasta file
      igi <- match(gi.taxid$V1[jtax], as.numeric(gi))
      # get the amino acid compositions of the sequences
      # use suppressMessages to hide mesages about unrecognized amino acid codes
      # set id="refseq" to avoid parsing IDs from fasta header lines
      aa <- suppressMessages(read.fasta(file="", i=ihead[igi], lines=fa, ihead=ihead, id="refseq"))
      # the number of amino acids read from file
      naa <- sum(aa[, 6:25])
      # prepare the dataframe
      # average amino acid composition
      aa[1, 6:25] <- colSums(aa[, 6:25])/nrow(aa)
      aa <- aa[1, ]
      aa$organism <- utax[j]
      aa$chains <- 1
      aa$abbrv <- naa
      # what is the name of this organism: look first at multispecies WP, then in gi header
      jtax.WP <- match(utax[j], gi.taxid.WP$V2)
      if(!is.na(jtax.WP)) orgname <- paste("[", gi.taxid.WP$V3[jtax.WP], "]", sep="")
      else {
        orgname <- s2c(fa[ihead[igi[1]]], sep=" [", keep.sep=TRUE)[2]
        orgname <- substr(orgname, 2, nchar(orgname))
      }
      cat(paste(orgname, "\n"))
      # append columns for numbers of sequences, organism name
      aa <- cbind(aa, data.frame(nseq=nseq, orgname=orgname))
      # keep track of source file name, and number of sequences and amino acids
      fname <- sub("microbial.nonredundant_protein", "nonredundant", sub(".protein.faa.gz", "", files[i]))
      aa$ref <- paste(fname, "(", nseq, ",", naa, ")", sep="")
      if(j==1) aa.out <- aa else aa.out <- rbind(aa.out, aa)
    }
    if(is.null(out)) out <- aa.out else out <- rbind(out, aa.out)
  }

  # now we have to combine taxids that showed up more than once (in different files)
  duptax <- unique(out$organism[duplicated(out$organism)])
  if(length(duptax) > 0) {
    for(i in 1:length(duptax)) {
      id <- which(out$organism==duptax[i])
      # add them together, weighted by numbers of sequences
      aa <- colSums(out$nseq[id]*out[id, 6:25]) / sum(out$nseq[id])
      # keep the result in the first row
      out[id[1], 6:25] <- aa
      # total number of sequences found in the files
      out$nseq[id[1]] <- sum(out$nseq[id])
      # total number of amino acids found in the sequence files
      out$abbrv[id[1]] <- sum(out$abbrv[id])
      # write down the file names 
      out$ref[id[1]] <- c2s(out$ref[id], sep=";")
      # drop the rows after the first
      out <- out[-id[2:length(id)], ]
    }
  }
  # append the organism name to the numbers of sequences and amino acids read
  for(i in 1:nrow(out)) {
    out$ref[i] <- paste(out$ref[i], out$orgname[i], sep="")
  }
  # we can drop the extra columns now
  out <- out[, 1:25]
  # and do a little rounding
  out[, 6:25] <- round(out[, 6:25], 3)
  write.csv(out, "protein_refseq_complete.csv", row.names=FALSE)
  return(out)
}
