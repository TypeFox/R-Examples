build.isotopeframe <-
function(pred.file, prey.file, prey.header = TRUE, pred.header = FALSE, sep="\t") {
  ## given a local filename leading to a tab-delimited file of isotope values, read in the
  ## values from the file and store them in a dataframe
  ## preyfiles is a vector of filenames, each containing a list of prey isotope measurements
  ## predatorfiles is a single filename, containing a list of predator isotope measurements
  ## sep is the delimiter used in the files. We assume tab delimited files
  prey.frame = read.table(prey.file, header=prey.header, sep=sep)
  pred.frame = read.table(pred.file, header=pred.header, sep=sep)
  list(prey.frame=prey.frame, pred.frame=pred.frame)
}
