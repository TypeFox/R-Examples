
readFreesurferAsciiHeader <- function(fileName) {
  ninfo <- as.integer(strsplit(readLines(fileName, n=2)[2], " ")[[1]])
  list(vertices=ninfo[1], faces=ninfo[2], label=stripExtension(FREESURFER_ASCII_SURFACE_DSET, basename(fileName)), embedDimension=3, headerFile=fileName, dataFile=fileName)
}

#' @importFrom readr read_table
readFreesurferAsciiGeometry<- function(fileName) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Pkg needed for this function to work. Please install it.",
         call. = FALSE)
  }
  ninfo <- as.integer(strsplit(readLines(fileName, n=2)[2], " ")[[1]])
  asctab <- read_table(fileName, skip=2)
  
  vertices <- as.matrix(asctab[1:ninfo[1],1:3])
  nodes <- as.matrix(asctab[(ninfo[1]+1):nrow(asctab),1:3])
  
  list(mesh=rgl::tmesh3d(vertices, nodes), headerFile=fileName, dataFile=fileName)
  
}

