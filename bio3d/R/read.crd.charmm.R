"read.crd.charmm" <-
function(file, ext=TRUE, verbose = TRUE, ...) {

  split.string <- function(x) {
    x <- substring(x, first, last)
    x[nchar(x) == 0] <- as.character(NA)
    x
  }
  
  trim <- function(s) {
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s == "")] <- NA
    s
  }

  if(ext)
    atom.format <- c(10, 10, 8, 8,
                     -4, 20,20,20, -1,
                     8, -1, 8, 20)
  else
    atom.format <- c(5, 5, 4, 5,
                     -1, 10,10,10, -1,
                     4, -1, 4, 10)
  
  atom.names  <- c("eleno", "resno", "resid", "elety",
                   "blank", "x", "y", "z",  "blank",
                   "segid", "blank", "resno2", "b")
  
  widths <- abs(atom.format)
  drop.ind <- (atom.format < 0)
  
  st <- c(1, 1 + cumsum(widths))
  first <- st[-length(st)][!drop.ind]
  last <- cumsum(widths)[!drop.ind]
  
  raw.lines <- readLines(file)

  
  head.ind <- which(substr(raw.lines,1,1)=="*")
  head.ind <- c(head.ind, (head.ind[length(head.ind)]+1) )

  if(length(head.ind)>0) {
    raw.lines <- raw.lines[-head.ind]
    if(verbose)
      cat(raw.lines[head.ind],sep="\n")
  }
  atom <- matrix(trim(sapply(raw.lines, split.string)),
                 byrow = TRUE,
                 ncol = length(atom.format[!drop.ind]),
                 dimnames = list(NULL,
                   atom.names[!drop.ind]) )
  
  output <- list(atom = atom,
                 xyz = as.numeric(t(atom[, c("x", "y", "z")])),
                 calpha = as.logical(atom[, "elety"] == "CA"))
  class(output) <- c("charmm", "crd")
  return(output)
}

