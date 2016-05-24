################################################################################
# File reading functions: read comparison
################################################################################
# read comparison file. Use source=tab or blast to specify type
read_comparison_from_blast <- function(file, sort_by="per_id",
                                       filt_high_evalue=NULL,
                                       filt_low_per_id=NULL,
                                       filt_length=NULL,
                                       color_scheme=NULL, ...){
  # read table anyway
  table <- read.table(file, as.is=TRUE, header=FALSE, sep="\t", quote="")
  # from blast output
  col_names <- c("name1", "name2", "per_id", "aln_len", "mism", "gaps",
                 "start1", "end1", "start2", "end2", "e_value", "bit_score")
  names(table) <- col_names
  # sort?
  if (!is.null(sort_by)){
    if (!sort_by %in% col_names) stop("Argument sort_by not in col names")
    if (!is.numeric(table[[sort_by]]))
      stop("Argument sort_by must designate a numeric column")
    decr <-
      if (sort_by %in% c("per_id", "aln_len", "bit_score")) FALSE else TRUE
    # reorder from weakest to strongest
    table <- table[order(table[[sort_by]], decreasing=decr),]
  }
  # filter evalue
  if (!is.null(filt_high_evalue)){
    if (!is.numeric(filt_high_evalue)) stop("filt_high_evalue must be numeric")
    table <- table[table$e_value <= filt_high_evalue,]
  }
  # filter per id
  if (!is.null(filt_low_per_id)){
    if (!is.numeric(filt_low_per_id)) stop("filt_low_per_id must be numeric")
    table <- table[table$per_id >= filt_low_per_id,]
  }
  # reorder columns
  table <- table[,c(col_names[c(7:10,1:6,11:12)])]
  
  # send further
  comparison <- .read_comparison(table, filt_length=filt_length, ...)
  # apply color scheme if requested
  if (!is.null(color_scheme)){
    comparison$col <- apply_color_scheme(comparison$per_id,
                                         direction=comparison$direction,
                                         color_scheme=color_scheme)
  }
  comparison
}
read_comparison_from_tab <- function(file, header=TRUE, ...){
  col_names <-  c("start1", "end1", "start2", "end2", "col")
  # from tabular data
  table <- read.table(file, as.is=TRUE, header=header, sep="\t", quote="")
  if (ncol(table) < 4) {
    stop("Insufficent number of columns in table")
  } else if (!header){
    if (ncol(table) == 4) {
      names(table) <- col_names[1:4]
    } else {
      names(table)[1:length(col_names)] <- col_names[1:length(col_names)]
    }
  }
  .read_comparison(table, ...)
}
.read_comparison <- function(table, filt_length=NULL, reverse=0,
                             color_scheme=NULL, ...){
  # check arguments
  if (ncol(table) < 4) stop("Insufficent number of columns in table")
  if (!all(c("start1", "end1", "start2", "end2") %in% names(table)))
    stop("Table must contain columns start1, end1, start2, end2")
  
  # filter length
  if (!is.null(filt_length)){
    if (!is.numeric(filt_length)) stop("filt_length must be numeric")
    av_len <- apply(cbind(abs(table$end1-table$start1),
                          abs(table$end2-table$start2)), 1, mean)
    table <- table[av_len >= filt_length,]
  }
  # reverse
  if (is.numeric(reverse) && reverse > 0){
    table <- reverse.comparison(table, side=reverse)
  }
  # make comparison object
  as.comparison(table)
}
