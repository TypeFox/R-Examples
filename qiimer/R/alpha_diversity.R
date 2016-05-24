#' Read a collated alpha diversity table from QIIME.
#'
#' @param filepath Path to collated alpha diversity file.
#' @return A data frame with columns `SampleID`, `sequences_per_sample`,
#'   `iteration`, and `diversity`.
#' @export
read_qiime_rarefaction <- function(filepath) {
  r <- read.table(
    filepath,
    sep="\t",
    quote="",
    na.strings=c('NA', 'n/a'),
    comment.char="#",
    header=TRUE)
  # Don't need to keep column of file names
  r <- r[2:length(colnames(r))]
  # Standardize column of rarefaction depth
  colnames(r)[1] <- "sequences_per_sample"
  sample_ids <- colnames(r)[3:length(colnames(r))]
  r <- reshape(
    r, 
    direction="long", 
    varying=sample_ids, 
    v.names="diversity", 
    timevar="SampleID",
    times=sample_ids)
  rownames(r) <- 1:nrow(r)
  r$SampleID <- factor(r$SampleID)
  r[,c("SampleID", "sequences_per_sample", "iteration", "diversity")]
}

#' Compute summary statistics for collated alpha diversity tables.
#'
#' @param rarefaction_table A collated alpha diversity data frame.
#' @return A data frame of summary statistics, with columns `SampleID`, 
#'   `sequences_per_sample`, `diversity.mean`, and `diversity.sd`.
#' @export
#' @examples 
#' data(relmbeta_alpha)
#' head(rarefaction_stats(relmbeta_alpha))
rarefaction_stats <- function(rarefaction_table) {
  # Aggregating with 2 functions makes the dimensions wonky.  We want
  # a flat table, so cast result to list and then back to data frame.
  as.data.frame(as.list(aggregate(
    diversity ~ SampleID + sequences_per_sample,
    data=rarefaction_table,
    FUN=function(x) c(mean=mean(x), sd=sd(x)))))
}
