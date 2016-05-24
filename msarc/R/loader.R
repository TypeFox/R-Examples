# Read nested table, get colnames from 1 header line, skipping sub-lines
# (blank in first column).  If the file doesn't have at least one data line,
# the function will fail.
read.nested <- function(filename) {
  conn <- file(filename,open="r")
  lines = readLines(conn,n=-1,warn=F)
  close(conn)

  # make the data frame... colnames from first line, minus trailing tabs
  lines <- sub("\t*$","",lines)
  fields <- strsplit(lines,"\t")
  names <- gsub("^\\s+|\\s+$","",fields[[1]])
  first <- fields[[2]]
  frame <- data.frame(t(first),stringsAsFactors=F)
  colnames(frame) <- names

  # read the rest, looking for "main" lines
  for (i in 3:length(lines)) {
    if (grepl("^[[:space:]]*$",lines[i])) {
      # skip blank line
    } else if (fields[[i]][1] != "") {
      frame <- rbind(frame,fields[[i]])
    }
  }
  return(frame)
}

# extract symbol names from ms file descriptions.  Looks for a pattern of
# the form "- [XXX_YYYY]" where "XXX" is the symbol and "YYYY" is the species
# (though it doesn't check or use the species in any way)
extract.symbols <- function(descriptions) {
  return(sub("^.*-\\s\\[(\\w+)_[A-Z0-9]+\\]\"?","\\1",descriptions,perl=T))
}

# Load a table from a mass spec experiment.  If there is no column named
# "Heavy/Light", onlyDiff is ignored.  If there is such a column (i.e. this
# is a Silac experiment), then onlyDiff=T means only load lines which have
# both heavy and light versions (i.e. those we can use in a Silac analysis),
# otherwise (onlyDiff=F) loads every line, not just the lines which represent
# things identified in both heavy and light cases.  Mostly we only want the
# cases where we see both, because that's what Silac is for.  We expect to
# see columns labelled "Accession", "Description" and "Score".
msarc.loadMS <- function(file,onlyDiff=T,config=list()) {
  conf <- config.defaults
  conf[names(config)] = config
  if (length(grep("xlsx",file))>0) {
    tbl <- readWorksheetFromFile(file,sheet=1)
    tbl <- tbl[!is.na(tbl[[conf$accessionCol]]),]
  } else {
    tbl <- read.nested(file)
  }
  if (conf$heavyCol %in% colnames(tbl) && onlyDiff) {
    tbl_diff <- tbl[[conf$heavyCol]] != ""
    tbl <- tbl[tbl_diff,]
  }
  frame <- data.frame(tbl[[conf$accessionCol]],
                      extract.symbols(tbl[[conf$descriptionCol]]),
                      tbl[[conf$scoreCol]],
                      stringsAsFactors=F)
  obj = msarc(frame)
  obj$filename <- basename(file)
  return(obj)
}
