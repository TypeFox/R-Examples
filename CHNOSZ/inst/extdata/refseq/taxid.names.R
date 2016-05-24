# generate a table with names of
# of species,genus,family,order,class,phylum,superkingdom
# for each of the microbial taxa in RefSeq database

# uses functions in CHNOSZ to process taxonomy files
require(CHNOSZ)

# change this to the location where names.dmp and nodes.dmp are located
taxdir <- "./taxdump"

# get the taxids from protein_refseq.csv
pr <- read.csv("protein_refseq.csv")
taxid <- pr$organism

# read in the names and nodes
if(!exists("taxnames")) {
  cat("reading names...\n")
  taxnames <- getnames(taxdir)
}
if(!exists("taxnodes")) {
  cat("reading nodes...\n")
  taxnodes <- getnodes(taxdir)
}

# what ranks we want to get
ranks <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")

# start with an empty list
out <- rep(list(character()), length(ranks))
names(out) <- c(ranks)

# loop over taxids
ii <- seq_along(taxid)
for(i in ii) {
  # test if data are available for this taxid
  if(!taxid[i] %in% taxnames$id) {
    for(j in 1:length(ranks)) out[[j]] <- c(out[[j]],NA)
    next
  }
  # get taxids of all parents
  pids <- allparents(taxid[i], taxdir, nodes=taxnodes)
  # get ranks of all parents
  pranks <- getrank(pids, taxdir, nodes=taxnodes)
  # find which parents are in the required ranks
  ip <- match(ranks, pranks)
  # get names of these parents
  pnames <- sciname(pids[ip], taxdir, names=taxnames)
  # add results to output list
  for(j in 1:length(ranks)) out[[j]] <- c(out[[j]], pnames[j])
  # report progress
  if(i %% 50 == 0) cat(paste(i, ".. "))
}
# finish progress report
cat("done!\n")

# write results to a file
out <- as.data.frame(out)
out <- cbind(data.frame(taxid=taxid[ii], out))
write.csv(out, "taxid_names.csv", row.names=FALSE, quote=FALSE)
