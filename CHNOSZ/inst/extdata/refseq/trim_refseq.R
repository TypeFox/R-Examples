# trim the protein_refseq.csv,
# keeping only organisms used in an example in ?protein.info
# (to keep the file size down for CHNOSZ package)

# the original file
pr <- read.csv("protein_refseq_complete.csv")
# the terms we'll keep
terms <- c("Natr", "Halo", "Rhodo", "Acido", "Methylo",
  "Nitro", "Desulfo", "Chloro", "Geo", "Methano",
  "Thermo", "Pyro", "Sulfo", "Buchner")
# identify rows to keep
iterm <- unique(unlist(sapply(terms, grep, pr$ref)))
# extract data frame, write new file
pr <- pr[iterm, ]
write.csv(pr, "protein_refseq.csv", row.names=FALSE)
