#' Parse an OTU mapping file from QIIME.
#'
#' @param filepath Path to the OTU mapping file.
#' @param prefix OTU names will be prefixed with this value.
#' @return A list of sequence identifiers for each OTU.
#' @export
read_qiime_otu_mapping <- function(filepath, prefix="") {
  # Based on http://stackoverflow.com/questions/6602881
  # Read in the data to character vector
  x <- scan(filepath, what="", sep="\n", quiet=T)
  # Separate elements by one or more whitepace
  y <- strsplit(x, "[[:space:]]+")
  # Extract the first vector element and set it as the list element name
  names(y) <- sapply(y, `[[`, 1)
  # names(y) <- sapply(y, function(x) x[[1]]) # same as above
  names(y) <- paste(prefix, names(y), sep="")
  # Remove the first vector element from each list element
  y <- lapply(y, `[`, -1)
  #y <- lapply(y, function(x) x[-1]) # same as above
  y
}

#' Create an OTU table.
#' 
#' @param otus A list of QIIME-format sequence identifiers for each OTU.
#' @param sample_ids An optional vector of sample IDs to include in the result.
#' @return A matrix of OTU counts by sample.
#' @export
make_otu_table <- function(otus, sample_ids=NULL) {
  sample_vec <- factor(sub("_\\d+$", "", unlist(otus), perl=T))
  otu_vec <- factor(rep(names(otus), sapply(otus, length)))
  if (!is.null(sample_ids)) {
    extra_levels <- setdiff(levels(sample_vec), sample_ids)
    if (length(extra_levels) > 0) {
      extra_idx <- which(sample_vec %in% extra_levels)
      sample_vec <- factor(sample_vec[-extra_idx], levels=sample_ids)
      otu_vec <- factor(otu_vec[-extra_idx])
    } else {
      sample_vec <- factor(sample_vec, levels=sample_ids)
    }
  }
  table(otu_vec, sample_vec, dnn=c("OTU", "SampleID"))
}
