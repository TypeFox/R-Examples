# Internal functions -----------------------------------------------------------

# Compute data for a specified gene.
# If specified, a percentage of gene values are deleted. By
# default, all the data are counted.
#' @import data.table
#' @importFrom utils count.fields
.gene_coverage_from_reads <-  function(geneSampleTable,
                                       percentage, sample = FALSE,
                                       introns, sf) {
  .shuffle <- function(file, percentage) {
    n <- length(count.fields(file)) * (1 - percentage)
    data <- tryCatch(
      data.table::fread(file),
      error = function(e) stop("From .gene_coverage_from_reads function: ", e)
    )
    data[sample(nrow(data), n)]
  }


  if (sample) {
    reads <- lapply(geneSampleTable$files,
                    .shuffle, percentage = percentage)
  }else{
    reads <- geneSampleTable$files
  }

  gene_data <- .get_data_from_reads(reads, geneSampleTable$gene_information,
                               geneSampleTable$types, sf)

  types <- geneSampleTable$types
  type1_u <- unique(types)[1]
  type1 <- paste0(l <- types[which(types == type1_u)],
                  "_rep", 1:length(l))
  type2_u <- unique(types)[2]
  type2 <- paste0(types[l <- which(types == type2_u)],
                  "_rep", 1:length(l))

  gene_data$mean_type1 <- rowMeans(gene_data[, type1, with = F])
  gene_data$mean_type2 <- rowMeans(gene_data[, type2, with = F])

  ## remove introns regions
  if (!is.null(introns)) {
    clean_data <- .remove_introns(data = gene_data,
                                  name = geneSampleTable$gene_information$name,
                                  introns = introns)
  }else{
    original_pos <- gene_data$pos
    gene_data$pos <- seq_len(length(gene_data$chrom))
    clean_data <- list(original_pos = original_pos,
                       data = gene_data,
                       introns = NULL)
  }

  gene_data <- as.list(clean_data$data)

  if (gene_data$strand[1] == "-") {
    gene_data <- lapply(gene_data[c("pos", "mean_type1", "mean_type2")], rev)
  }else{
    gene_data <- gene_data[c("pos", "mean_type1", "mean_type2")]
  }
  gene_data$original_pos <- clean_data$original_pos

  gene_data$gene_information <- geneSampleTable$gene_information
  gene_data$introns <- clean_data$introns
  gene_data$samples = list(type1 = type1, type2 = type2)

  class(gene_data) <- "geneCoverage"
  gene_data
}


# Calculates the f score using the bootstrap approach.
# @param boostrapnumber The number of time to reproduce the experiment
# @param reads The reads
# @param introns
# @param sf
# @param percentage
# @return A vector contening the distance of each point to the diagonale
.fscore_bootstrap <- function(boostrapnumber, geneSampleTable,
                              introns, sf, percentage){

  data <- .gene_coverage_from_reads(
    geneSampleTable, percentage,
    sample = TRUE,
    introns, sf
  )

  res <- .perpendicular_distance(
    xValues = data$pos,
    yValues = .difference_of_cdf(data[2:3])
  )
  rm(data); gc()
  return(res)
}

## Methods section ----------------

## Calulate multiple ranges of the data.
.calculate_intervals <- function(vector){
  isNext <- function(position, vector){
    ifelse( position == 1, TRUE, ifelse(vector[position] == vector[position - 1] +1, TRUE, FALSE))
  }
  l <- length(vector)
  intervalPoints <- which( sapply(c(1:l),isNext, vector= vector) == FALSE)
  nbpoints <- length(intervalPoints)+1
  intervalPoints <- c(1, intervalPoints-1, intervalPoints, l)

  res <- NULL
  for(i in 1:nbpoints) {
    res <- rbind(res, c(start = vector[intervalPoints[2*i-1]], end = vector[intervalPoints[2*i]]))
  }
  rownames(res) <- paste0("zone", 1:nbpoints)
  res
}

#calculates positions for which the observed f score is in the distribution of simulated f max.
.methodA <- function(observed_distances, q, selected){
  if(selected){
    tss_positions <- which(observed_distances > q[1] & observed_distances < q[2])
    .calculate_intervals(tss_positions)
  }
}

#calculates positions for which the mean simulated f score is in the distribution of simulated f max.
.methodB <- function(simulated_scores, q, selected){
  if(selected){
    mean <- apply(do.call(cbind, simulated_scores), 1, mean)
    tss_positions <- which(mean > q[1] & mean < q[2])
    .calculate_intervals(tss_positions)
  }
}

#calculates positions for which the simulated f score is in the distribution of simulated f max.
.methodC <- function(simulated_scores, q, selected){
  if(selected){
    tss_positions <- unlist(lapply(simulated_scores, function(value){
      which(value > q[1] & value < q[2])
    }))

    unique_positions <- unique(sort(tss_positions))
    .calculate_intervals(unique_positions)
  }
}

#calculates positions for which the simulated f score is in the distribution of simulated f max and
#calculates the gaussian mean and standrad deviation
#' @import mclust
.methodC_gaussian <- function(simulated_scores, q, selected){

  if(selected){
    tss_positions <- unlist(lapply(simulated_scores, function(value){
      which(value > q[1] & value < q[2])
    }))

    improoved_gaussians <- function(g, last_clusters){
      if(g > 10){
        return (last_clusters)
      }
      suppressWarnings(clusters <- mclust::Mclust(tss_positions , G = 1:g))
      improvement <- abs(clusters$bic - last_clusters$bic) /
        abs(last_clusters$bic) * 100
      ifelse(improvement > 1,
             return(improoved_gaussians(g+1, clusters)),
             return(last_clusters))
    }

    gaussians_g1 <- mclust::Mclust(tss_positions, G = 1)
    gaussians <- improoved_gaussians (g = 2, last_clusters = gaussians_g1)
    mean_values <- gaussians$parameters$mean
    sd_values <- sqrt(gaussians$parameters$variance$sigmasq)
    list(mean = round(mean_values), sd = round(sd_values))
  }
}

#Calculates positions for each simulated f max
.methodD <- function(simulated_scores, selected){
  if(selected){
    tss_positions <- unique(unlist(lapply(simulated_scores, which.max)))
    .calculate_intervals(sort(tss_positions))
  }
}
#

# BAM to BED methods -------------------

#' Convert BAM files into BED files containing the reads overlaping the specified annotations.
#'
#'
#' @param bamfile a character vector indicating the BAM file name.
#'        Note: The bamfile must be sorted by coordinates.
#' @param annotations an object of type \code{\link{annotationsSet}}
#'        containing information on one genes.
#' @param paired_end logical indicating whether the \code{bamfile}
#'        contains paired-end data.
#' @param as_fragments logical indicating if paired-end data must be paired
#'        and merged to form fragments.
#' @param outfile a character vector indicating the output file name.
#'        If not provided, the result will be internalized in R.
#' @param flanking_region Number of bases before and after the gene ORF should be included.
#'
#' @return An object of type \code{\link[data.table]{data.table}} with 9 columns.
#'  \tabular{ll}{
#'    chromosome        \tab Read chromosome.\cr
#'    start             \tab Read starting position.\cr
#'    end               \tab Read ending position.\cr
#'    name              \tab Query template name .\cr
#'    score             \tab Mapping quality.\cr
#'    ChromosomeNext    \tab Gene strand.\cr
#'    startNext         \tab Position of the mate/next read.\cr
#'    tlen              \tab observed read length.\cr
#'    flag              \tab bitwise flag.
#'  }
#'
#' @export
#' @examples
#' bamfile <- system.file("extdata", "wt_rep1.bam", package = "yCrypticRNAs")
#' data(annotations)
#' bam_to_reads(bamfile, annotations)
bam_to_reads <- function(bamfile, annotations, paired_end = TRUE,
                         as_fragments = TRUE, outfile = NULL, flanking_region = 0) {

  if(!paired_end & as_fragments){
    stop("You can't make fragments for data that are not paired-end")
  }

  if(!file.exists(bamfile)){
    stop("The file provided doesn't exists.")
  }

  index <- paste0(bamfile, ".bai")
  unlink_index <- .index_bam(bamfile, index)

  if(!is.null(outfile))
    if(file.exists(outfile))
      unlink(outfile)

  genes <- .get_genes(annotations)
  res <- do.call(rbind, lapply(genes, function(gene) {
    annotation <- .get_annotation(annotations, gene)
    .bam_to_reads_by_gene(bamfile, annotation, paired_end,
                          as_fragments, outfile, flanking_region)
  }))

  if(unlink_index)
    unlink(index)
  if(is.null(outfile))
    return (res)
}

.index_bam <- function(infile, outfile){
  unlink_index <- FALSE
  if(!file.exists(outfile)){
    unlink_index <- TRUE
    suppressWarnings(tryCatch(
      Rsamtools::indexBam(infile),
      error = function(e){
        tryCatch(
          Rsamtools::sortBam(infile, substr_left(infile, -4)),
          error = function(e2){
            stop(e2)
          }
        )
        Rsamtools::indexBam(infile, outfile)
      }
    ))
  }
  return(unlink_index)
}

.bam_to_coverage <- function(bamfile, annotations, sf, paired_end = TRUE,
                         as_fragments = TRUE, outfile = NULL) {

  index <- paste0(bamfile, ".bai")
  unlink_index <- .index_bam(bamfile, index)

  if(!is.null(outfile))
    if(file.exists(outfile))
      unlink(outfile)

  genes <- .get_genes(annotations)

  progress <- FALSE
  if(length(genes) > 1){
    print(paste(bamfile, ":"))
    pb <- txtProgressBar(min = 0, max = length(genes), label = bamfile, style = 3)
    progress <- TRUE
  }
  res <- do.call(rbind, lapply(genes, function(gene) {
    annotation <- .get_annotation(annotations, gene)

    if(progress)
      setTxtProgressBar(pb, which(genes == gene))

    reads <- .bam_to_reads_by_gene(bamfile, annotation, paired_end, as_fragments, outfile, flanking_region = 100)
    coverage(reads, annotation, sf)
  }))
  if(progress)
    close(pb)

  if(unlink_index)
    unlink(index)
  if(is.null(outfile))
    return (res)
}

.bam_to_reads_by_gene <- function(bamfile, annotation, paired_end,
                                  as_fragments, outfile, flanking_region) {


  gene <- IRanges::RangesList(IRanges::IRanges(start = annotation$start - flanking_region,
                                               end = annotation$end + flanking_region))
  names(gene) <- annotation$chromosome
  gene_bam_file <- tempfile()

  Rsamtools::filterBam(bamfile, gene_bam_file,
                       param = Rsamtools::ScanBamParam(which = gene))

  if (annotation$strand == "-")
    reads <- .select_neg_strand(gene_bam_file, paired_end, as_fragments)
  else
    reads <- .select_pos_strand(gene_bam_file, paired_end, as_fragments)

  unlink(gene_bam_file)
  if (is.null(outfile))
    return(reads)
  .writeFile(reads, outfile, append = T)
}

.select_flag <- function(bamfile, flag){

  what <- c("qname", "flag", "rname", "pos", "qwidth", "mapq", "mate_status", "strand")
  sam_data <- Rsamtools::scanBam(
    bamfile,
    param = Rsamtools::ScanBamParam(flag = flag, what = what)
  )
  .samtobed(sam_data)
}

.samtobed <- function(sam_data){
  sam_data <- sam_data[[1]]
  bed_data <- data.table::data.table(
    chromosome = sam_data$rname,
    start = sam_data$pos,
    end = sam_data$pos + sam_data$qwidth,
    name = sam_data$qname,
    mapq = sam_data$mapq
  )
  bed_data <- bed_data[ bed_data$mapq == 50]
  data.table::setorder(bed_data, name)
}

.as_fragments <- function(first, second){

  first_no_mate <- which(!is.element(first$name, second$name))
  second_no_mate <- which(!is.element(second$name, first$name))

  if(length(first_no_mate) > 0)
    first <- first[ - first_no_mate]
  if(length(second_no_mate) > 0)
    second <- second[ - second_no_mate]

  first[, end := second$end]
  return(first)
}

.select_neg_strand <- function(bamfile, paired_end, as_fragments){
  if(paired_end){

    neg99_flag <- Rsamtools::scanBamFlag(
      isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE, isMinusStrand = FALSE,
      isMateMinusStrand = TRUE, isFirstMateRead = TRUE,
      isSecondMateRead = FALSE, isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE, isDuplicate = NA
    )

    neg147_flag <- Rsamtools::scanBamFlag(
      isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE, isMinusStrand = TRUE,
      isMateMinusStrand = FALSE, isFirstMateRead = FALSE,
      isSecondMateRead = TRUE, isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE, isDuplicate = NA
    )
    neg99 <- .select_flag(bamfile, neg99_flag)
    neg147 <- .select_flag(bamfile, neg147_flag)

    if (is.null(neg99) || is.null(neg147))
      return (NULL)
    if (as_fragments){
      res <- .as_fragments(neg99, neg147)
      res [, strand := "-"]
      return (res)
    }

    res <- merge(neg99, neg147, by = names(neg99), all = T)
    res [, strand := "-"]
    return (res)
  }
  flag <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE, isMinusStrand = TRUE,
    isSecondaryAlignment = FALSE, isDuplicate = FALSE,
    isNotPassingQualityControls = FALSE
  )
  res <- .select_flag(bamfile, flag)
  res [, strand := "-"]
  return(res)
}

.select_pos_strand <- function(bamfile, paired_end, as_fragments){
  if(paired_end){
    pos83_flag <- Rsamtools::scanBamFlag(
      isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE, isMinusStrand = TRUE,
      isMateMinusStrand = FALSE, isFirstMateRead = TRUE,
      isSecondMateRead = FALSE, isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE, isDuplicate = NA
    )
    pos163_flag <- Rsamtools::scanBamFlag(
      isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE, isMinusStrand = FALSE,
      isMateMinusStrand = TRUE, isFirstMateRead = FALSE,
      isSecondMateRead = TRUE, isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE, isDuplicate = NA
    )

    pos83 <- .select_flag(bamfile, pos83_flag)
    pos163 <- .select_flag(bamfile, pos163_flag)

    if(is.null(pos83) || is.null(pos163))
      return ( NULL)
    if(as_fragments){
      res <- .as_fragments(pos163, pos83)
      res[, strand := "+"]
      return (res)
    }

    res <- merge(pos83, pos163, by = names(pos83), all = T)
    res[, strand := "+"]
    return (res)
  }
  flag <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE, isMinusStrand = TRUE,
    isSecondaryAlignment = FALSE, isDuplicate = FALSE,
    isNotPassingQualityControls = FALSE
  )
  res <- .select_flag(bamfile, flag)
  res[, strand := "+"]
  return(res)
}

# cTSS function --------------------

#' @title Cryptic transcription start sites (cTSS)
#'
#' @description For genes expresing a cryptic transcript, the method allows the
#'   identification of cTSS. Five differents methods are available.
#'
#' @param bamfiles a vector of characters indicating the BAM file paths.
#' @param name The name of the gene for which to calculate the cryptic
#'   initiation site(s).
#' @param annotations An object of type \code{\link{annotationsSet}}
#'        containing information on genes.
#' @param introns an objet of type \code{\link{annotationsSet}}
#'        containing the annotations of the intronic regions.
#'        Note: The introns must have same name as the gene they
#'        are associated with.
#' @param replicates The number of time to sample the data. Default = 200.
#' @param percentage Th fraction of data to be removed at each simulation.
#'        Default = 0.1.
#' @param sf A vector of the scaling factors to apply to each
#'        sample. Must be the same length as the \code{fragments_file}
#' @param paired_end logical indicating whether the \code{bamfile}
#'        contains paired_end data.
#' @param as_fragments logical indicating if paired_end data must paired
#'        and merged to form fragments.
#' @param method A character string or a vector specifying the method to use to
#'   calculate the cryptic initiations sites. Must be one of "methodC_gaussian"
#'   (default), "methodA", "methodB", "methodC" or "methodD".
#' @param types a vector of the same length as bamfiles indicating the type
#'        of data in each file. Example: WT vs MUT or untreated vs treated".
#'
#' @details
#' By definition, the observed f value for a gene is the perpendicular
#'  distance between the differential cumulative RNA-seq values (type1 - type2)
#'  and a diagonal linking the first and last data points. The simulated
#'  f max is the maximum f value for a gene after re-sampling the data.
#'
#' Method A identify a cryptic zone by calculating positions for which
#' the observed f value is in the distribution of simulated f max.
#'
#' Method B identify a cryptic zone by calculating positions for which the mean
#' simulated f value is in the distribution of simulated f max.
#'
#' Method C identify a cryptic zone by calculating positions for which the
#' simulated f value is in the distribution of simulated f max.
#'
#' Method D identify a cryptic zone by calculating the positions for each
#' simulated f max.
#'
#' Method C gaussian determine the mean and standard deviation of all the
#' positions for which the simulated f value is in the distribution of simulated
#' f max.
#'
#' @export
#' @return A list with the following components:
#'       \item{methodC_gaussian}{cTSS mean and sd values using the
#'                               method C (gaussian)}
#'        \item{methodA}{cryptic zones start and end
#'                        position using the method A}
#'        \item{methodB}{cryptic zones start and end
#'                        position using the method B}
#'        \item{methodC}{cryptic zones start and end
#'                        position using the method C}
#'        \item{methodD}{cryptic zones start and end
#'                        position using the method D}
#'        \item{gene_information}{An object of class \code{\link{annotationsSet}}
#'                              containing the information on the gene.}
#'
#' @importFrom stats quantile
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data("annotations")
#' samples <- c("wt_rep1", "wt_rep2", "mut_rep1", "mut_rep2")
#' bamfiles <- system.file("extdata", paste0(samples, ".bam"),
#'                         package = "yCrypticRNAs")
#' sf <- c(0.069872847, 0.081113079, 0.088520251, 0.069911116)
#' data(introns)
#' types = c("wt", "wt", "mut", "mut")
#'
#'
#' initiation_sites("YER109C", bamfiles, types, annotations,
#'                  introns, sf = sf, replicates = 5)
#' #initiation_sites("YER109C", bamfiles, types, annotations,
#'                 # introns, sf = sf, percentage = 0.2, method = "methodA")
initiation_sites <- function(name, bamfiles, types, annotations, introns = NULL,
                             sf = replicate(length(bamfiles), 1), replicates = 200, percentage = 0.1,
                             method = c("methodC_gaussian", "methodA",
                                        "methodB", "methodC", "methodD"),
                             paired_end = TRUE, as_fragments = TRUE){

  method <- match.arg(method, several.ok = T)

  if(!paired_end & as_fragments)
    stop("You can't make fragments for data that are not paired-end")

  if(length(bamfiles) != length(sf))
    stop ("You must provide one file per sample and give a scaling factor for each file")

  if(!is.element( "annotationsSet", class(annotations))){
    annotations <- as.annotationsSet(annotations)
  }

  if(!is.element("annotationsSet", class(introns))){
    introns <- as.annotationsSet(introns)
  }

  annotation <- .get_annotation(annotations, name)
  if(is.element(name, annotations$name)){
    gene_reads <- paste0(substr_left(bamfiles,-4), "_", name, ".txt")
    mapply(
      bam_to_reads,
      bamfile = bamfiles,
      annotation = list(annotation),
      paired_end = paired_end,
      as_fragments = as_fragments,
      outfile = gene_reads,
      flanking_region = 100
    )
    empty_file = is.element(0, base::file.size(gene_reads))
  }else{
    empty_file <- TRUE
  }

  if(! empty_file){
    geneSampleTable <- list(gene_information = annotation,
                            types = types,
                            files = gene_reads)

    data <- .gene_coverage_from_reads(geneSampleTable = geneSampleTable,
                                      introns = list(introns), percentage = percentage,
                                      sf = sf, sample = F)

    cdf <- .difference_of_cdf (data = data[2:3])
    observed_distances <- .perpendicular_distance(xValues = data$pos, yValues = cdf)

    simulated_scores <- parallel::mclapply(as.list(1:replicates), .fscore_bootstrap,
                                           geneSampleTable = geneSampleTable,
                                           percentage = percentage, introns = introns,
                                           sf = sf)
    simulated_fmax <- base::sapply(simulated_scores, max)
    q <- quantile(simulated_fmax, probs = c(2.5, 97.5)/100)
    tss_positions <- which(observed_distances > q[1] & observed_distances < q[2])
  }else{
    tss_positions <- NULL
  }

  if( length(tss_positions) == 0 ){
    if(is.element(name, annotations$name))
      lapply(gene_reads, unlink)

    result <- list( NA, NA, NA, NA, NA)
    names(result) <- c("methodC_gaussian", "methodA", "methodB", "methodC", "methodD")
    result$gene_information <- annotation
    return(result)
  }

  resA = resB = resC = resC_gaussian = resD = NULL

  all_methods <- c("methodC_gaussian", "methodA", "methodB", "methodC", "methodD")
  selected_methods <- is.element(all_methods, method)

  result <- list( .methodC_gaussian(simulated_scores = simulated_scores, q = q, selected_methods[1]),
                  .methodA(observed_distances = observed_distances, q = q, selected_methods[2]),
                  .methodB(simulated_scores = simulated_scores, q = q, selected_methods[3]),
                  .methodC(simulated_scores = simulated_scores, q = q, selected_methods[4]),
                  .methodD(simulated_scores = simulated_scores, selected_methods[5]))

  names(result) <- c("methodC_gaussian", "methodA", "methodB", "methodC", "methodD")
  result$gene_information <- annotation

  #Delete unneeded files
  lapply(gene_reads, unlink)
  result
}
