#' Read a QIIME sample mapping file.
#'
#' @param filepath Path to sample mapping file.  The file must conform to the
#'   QIIME standards, detailed at
#'   \url{http://qiime.org/documentation/file_formats.html}.
#' @return A data frame of sample information.  Because the SampleID column is
#'   so often used to extract data from distance matrices and OTU tables, it 
#'   is returned as a character vector.
#' @export
read_qiime_mapping_file <- function(filepath) {  
  sample_file <- file(filepath, 'rt')

  header <- readLines(sample_file, n=1)
  header <- sub("^#", "", header)
  cols <- unlist(strsplit(header, "\t"))

  sample_data <- read.table(
    sample_file,
    col.names=cols,
    sep="\t",
    quote="",
    row.names=NULL,
    na.strings=c("NA", "na", "Null", "null"))
  close(sample_file)

  # The SampleID column is often used to extract data from distance matrices 
  # and OTU tables.  Storing as a character vector facilitates this practice.
  sample_data$SampleID <- as.character(sample_data$SampleID)
  
  # Ideally, the Description column should contain a unique free text
  # description of each sample.  In this case, there is no advantage to using
  # a factor data type.
  if ('Description' %in% names(sample_data)) {
    sample_data$Description <- as.character(sample_data$Description)
  }

  sample_data
} 
