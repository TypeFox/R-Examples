require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
seqs <- read.ms(seqFile)
# Print the number of sequences in the MS object and whether
#   stored in C or R
print(seqs)
