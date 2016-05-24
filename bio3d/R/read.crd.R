"read.crd" <- function(file, ...) {
  ## from tools package:
  pos <- regexpr("\\.([[:alnum:]]+)$", file)
  ext <- ifelse(pos > -1L, substring(file, pos + 1L), "")

  if(ext %in% c("crd")) {
    class(file)=c("character", "charmm")
  }

  if(ext %in% c("inpcrd", "rst")) {
    class(file)=c("character", "amber")
  }

  UseMethod("read.crd", file)
}

