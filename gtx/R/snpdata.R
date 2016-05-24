summary.snpdata <- function(object, ...) {
  stopifnot(is.snpdata(object))
  cat("SNP information", paste(names(object$snpinfo), collapse = ","), "for", nrow(object$snpinfo), "SNPs\n")
  cat("Genotype data for", nrow(object$data), "subjects\n")
  extra <- setdiff(names(object$data), c(paste(object$snpinfo$snp, object$snpinfo$coded.allele, sep = "_"),
                                          paste(object$snpinfo$snp, object$snpinfo$coded.allele, sep = "_")))
  if (length(extra) > 0) cat("Additional data", paste(extra, collapse = ","), "for", nrow(object$data), "subjects\n")
}

is.snpdata <- function(object) {
  if (class(object) != "snpdata") return (FALSE)
  if (!is.list(object)) return (FALSE)
  if (!all(c("snpinfo", "data") %in% names(object))) return (FALSE)
  if (!is.data.frame(object$snpinfo)) return (FALSE)
  if (!all(c("snp", "coded.allele", "noncoded.allele") %in% names(object$snpinfo))) return (FALSE)
  if (!is.data.frame(object$data)) return (FALSE)
  if (!all(paste(object$snpinfo$snp, object$snpinfo$coded.allele, sep = "_") %in% names(object$data))) return (FALSE)
  return (TRUE)
}

as.snpdata <- function(object) {
  if (!is.list(object)) return (NA)
  if (!all(c("snpinfo", "data") %in% names(object))) return (NA)
  if (!is.data.frame(object$snpinfo)) return (NA)
  if (!all(c("snp", "coded.allele", "noncoded.allele") %in% names(object$snpinfo))) return (NA)
  if (!is.data.frame(object$data)) return (NA)
  if (!all(paste(object$snpinfo$snp, object$snpinfo$coded.allele, sep = "_") %in% names(object$data))) return (NA)
  class(object) <- "snpdata"
  return (object)
}

