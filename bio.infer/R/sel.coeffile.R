# 1.7.2009
# Select coefficient file from those available

sel.coeffile <- function(fstring = "coef") {

  temp <- data(package = "bio.infer")
  flist <- temp$results[,3]

  incvec <- substring(flist, 1, 4) == substring(fstring,1,4)

  flist <- flist[incvec]

  fsel <- tk_select.list(flist, title = "Available coefficient files")

  return(fsel)
}

