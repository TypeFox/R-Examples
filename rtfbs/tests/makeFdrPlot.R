require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
ms <- read.ms(seqFile);
pwmFile <- "pwm.meme"
unzip(exampleArchive, pwmFile)
# Read in Position Weight Matrix (PWM) from MEME file from
#  the examples into a Matrix object
pwm <- read.pwm(pwmFile)
# Build a 3rd order Markov Model to represent the sequences
#   in the MS object "ms".  The Model will be a list of
#   matrices  corrisponding in size to the order of the 
#   Markov Model
mm <- build.mm(ms, 3);
# Match the PWM against the sequences provided to find
#   possible transcription factor binding sites.  A 
#   Features object is returned, containing the location
#   of each possible binding site and an associated score.
#   Sites with a negative score are not returned unless 
#   we set threshold=-Inf as a parameter.
cs <- score.ms(ms, pwm, mm)
# Generate a sequence 1000 bases long using the supplied
#   Markov Model and random numbers
v <- simulate.ms(mm, 10000)
# Match the PWM against the sequences provided to find
#   possible transcription factor binding sites.  A 
#   Features object is returned, containing the location
#   of each possible binding site and an associated score.
#   Sites with a negative score are not returned unless 
#   we set threshold=-Inf as a parameter. Any identified
#   binding sites from simulated data are false positives
#   and used to calculate False Discovery Rate
xs <- score.ms(v, pwm, mm)
# Calculate the False Discovery Rate for each possible
#   binding site in the Features object CS.  Return
#   a mapping between each binding site score and the
#   associated FDR.
fdr <- calc.fdr(ms, cs, v, xs)
# Plot the False Discovery Rate v.s. score for one or
#   more groups.  To plot multiple FDR/Score mapping 
#   data frames on a single plot, simply supply a list
#   of FDR/Score data.frames 
makeFdrPlot(fdr)
