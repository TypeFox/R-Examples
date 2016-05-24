# Default configuration values
config.defaults = list("heavyCol" = "Heavy/Light",
                       "accessionCol" = "Accession",
                       "descriptionCol" = "Description",
                       "scoreCol" = "Score",
                       "backgroundSaturation" = 25,
                       "backgroundLightness" = 99,
                       "radius" = 1000)

# Create an msarc object from a 3-column data frame containing the UniProt ID,
# the gene symbol, and the score (in that order).
msarc <- function(df) {
  obj = list()
  obj$config <- config.defaults
  colnames(df) <- c('uniprot','symbol','score')
  df$uniprot <- as.character(df$uniprot)
  df$symbol <- as.character(df$symbol)
  obj$rawdata <- df
  obj$data <- obj$rawdata
  obj$filename <- "data frame"
  obj$control <- "unset"
  obj$go2uni <- "unset"
  obj$gotbl <- "unset"
  obj$candidates <- "unset"
  obj$counts <- "unset"
  obj$tree <- "unset"
  class(obj) <- "msarc"
  return(obj)
}

# Subtract one MS list (or at least a list of Uniprot IDs) from
# a MS experiment.
msarc.subtract <- function(source,control) {
  if (class(control) == 'msarc') {
    lose <- source$data$uniprot %in% control$data$uniprot
    source$control <- control$filename
  } else {
    lose <- source$data$uniprot %in% control
    source$control <- "vector"
  }
  source$data <- source$data[!lose,]
  return(source)
}
  
summary.msarc <- function(object,...) {
  print(object)
  print("Sample of GO categories:")
  head(object$gotbl,n=20)
}

print.msarc <- function(x,...) {
  if (x$control == "unset") {
    descr <- sprintf("msarc:  proteins=%d\n  source=%s\n",
                      length(x$data$uniprot),
                      x$filename)
  } else {
    descr <- sprintf("msarc:  proteins=%d (%d pre-subtraction)\n  source=%s\n  control=%s\n",
                     length(x$data$uniprot),
                     length(x$rawdata$uniprot),
                     x$filename,
                     x$control)
  }
  if (class(x$go2uni) == "list") {
    if (length(x$go2uni) == length(x$go2uniAll)) {
      descr <- sprintf("%s  GO terms: %d (%d direct, %d by anscestral links)\n",descr,length(x$go2uni),length(x$go2uniBase),length(x$go2uniAll)-length(x$go2uniBase))
    } else {
      descr <- sprintf("%s  GO terms: %d (pre-filtering: %d direct, %d by ancestral links)\n",descr,length(x$go2uni),length(x$go2uniBase),length(x$go2uniAll)-length(x$go2uniBase))
    }
  }
  cat(descr)
  invisible(descr)
}
