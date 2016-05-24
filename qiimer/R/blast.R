#' Parse tabular output from BLAST.
#'
#' @param filepath Path to tabular BLAST output file.
#' @return A data frame of BLAST results.
#' @export
read_blast_table <- function (filepath) {
  column_classes <- c(
    "character", "character", "numeric",
    "integer", "integer", "integer",
    "integer", "integer", "integer", "integer",
    "numeric", "numeric"
  )
  column_names <- c(
    "query_id", "subject_id", "pct_identity", 
    "alignment_len", "mismatches", "gap_openings", 
    "query_start", "query_end", "subject_start", "subject_end", 
    "e_value", "bit_score")
  read.delim(
    filepath, header=F, comment.char="#", 
    colClasses=column_classes, col.names=column_names
  )
}
