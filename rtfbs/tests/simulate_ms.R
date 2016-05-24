require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
ms <- read.ms(seqFile);
# Build a 3rd order Markov Model to represent the sequences
#   in the MS object "ms".  The Model will be a list of
#   matrices  corrisponding in size to the order of the 
#   Markov Model
mm <- build.mm(ms, 3);
# Generate a sequence 1000 bases long using the supplied
#   Markov Model and random numbers
v <- simulate.ms(mm, 1000);
