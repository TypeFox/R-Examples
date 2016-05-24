utils::globalVariables(names = c("chrom", "start", "end", "name", "strand",
                                 "pos", "score", "V8", "V6"), package = "yCrypticRNAs")

#' Create a CoverageDataSet
#'
#' Create a CoverageDataSet containing the depth and breadth of coverage of genes in \code{annotations}.
#'
#' @param bamfiles A vector of characters indicating the BAM files paths.
#' @param annotations An object of type \code{\link{annotationsSet}}
#'        containing information on genes.
#' @param types A vector of the same length as \code{bamfiles} indicating the type
#'        of data in each file. Example: WT vs MUT or untreated vs treated".
#' @param sf A vector of the same length as \code{bamfiles} indicating the
#'        scaling factor to apply on each file. Each coverage value is multiplied
#'        by this factor before being reported. Useful for normalizing coverage
#'        by, e.g., reads per million (RPM).
#' @param paired_end logical indicating whether the \code{bamfiles}
#'        contains paired-end data.
#' @param as_fragments logical indicating if paired-end data must paired
#'        and merged to form fragments.
#'
#' @return An object of class 'coverageDataSet' containing the coverage for each sample.
#'         The details of the output componets are as follow:
#'
#'    \item{data}{ A \code{\link[data.table]{data.table}} with the following components:
#'      \describe{
#'        \item{chrom}{Gene chromosome.}
#'        \item{start}{Gene transcription start site.}
#'        \item{end}{Gene transcription termination site.}
#'        \item{name}{Gene name.}
#'        \item{score}{Gene score.}
#'        \item{strand}{Gene strand.}
#'        \item{pos}{One based positions of the gene.}
#'        \item{...}{For each sample, report the depth at each position of the gene per sample.}
#'      }}
#'    \item{samples}{A list containg the samples' names that are in the coverageDataSet.
#'     \describe{
#'        \item{type1}{Type1 samples names (controls)}
#'        \item{type2}{Type2 samples names (experiment)}
#'      }}
#'
#'
#' @export
#' @examples
#' samples <- c("wt_rep1", "wt_rep2", "mut_rep1", "mut_rep2")
#' bamfiles <- system.file("extdata", paste0(samples, ".bam"),
#'                               package = "yCrypticRNAs")
#' data(annotations)
#' types <- c("wt", "wt", "mut", "mut")
#' scaling_factors <- c(0.069872847, 0.081113079, 0.088520251, 0.069911116)
#' rna_seq_signals <- coverageDataSet (bamfiles, annotations, types, scaling_factors)
#' rna_seq_signals
coverageDataSet <- function(bamfiles, annotations, types,
                            sf = c(1,1,1,1), paired_end = TRUE,
                            as_fragments = TRUE) {

  if(length(bamfiles) != length(types))
    stop ("You must provide one file per sample and specify the type of each sample")

  if(length(bamfiles) != length(sf))
    stop ("You must provide one file per sample and give a scaling factor for each file")

  lapply(bamfiles, function(file){
    if(!file.exists(file) || file.info(file)$size == 0)
      stop(paste("The file:", file, "doesn't exist or is empty") )
  })

  if(!paired_end & as_fragments){
    stop("You can't make fragments for data that are not paired-end")
  }

  if(length(unique(types)) != 2)
    stop("The program handle only two types of samples. Example: WT vs MUT or untreated vs treated")

  # bam to cov or fragments
  cov <- mapply(.bam_to_coverage, bamfile = bamfiles,
                annotations = list(annotations),
                sf = sf, paired_end, as_fragments, SIMPLIFY = F)

  data <- .get_data_from_coverage(cov, types, sf)
  remove(cov)
  invisible(gc())

  type1_u <- unique(types)[1]
  type1 <- paste0(l <- types[which(types == type1_u)],
                  "_rep", 1:length(l))
  type2_u <- unique(types)[2]
  type2 <- paste0(types[l <- which(types == type2_u)],
                  "_rep", 1:length(l))

  res <- list(data = data,
              samples = list(type1 = type1, type2 = type2))

  class(res) <- append("coverageDataSet", class(res))
  res
}

.get_data_from_reads <- function(reads, annotations, types, sf){

  # reads to coverage
  cov <-  mapply(
    coverage, reads = reads,
    annotation = list(annotations),
    sf = sf, SIMPLIFY = F
  )
  .get_data_from_coverage( cov, types, sf)
}

.get_data_from_coverage <- function(cov, types, sf){

  type1_u <- unique(types)[1]
  type1 <- paste0(l <- types[which(types == type1_u)],
                  "_rep", 1:length(l))
  type2_u <- unique(types)[2]
  type2 <- paste0(types[l <- which(types == type2_u)],
                  "_rep", 1:length(l))

  type1_data <- cov[which(types == type1_u)]
  type2_data <- cov[which(types == type2_u)]

  .get_score<- function(data){
    data[,V8]
  }

  data <- data.table::data.table(type1_data[[1]][, 1:7, with = F],
                                 sapply(type1_data, .get_score),
                                 sapply(type2_data, .get_score))

  data.table::setnames(data, c("chrom", "start", "end", "name", "score", "strand", "pos", type1, type2))
  data.table::setkey(data, "name")
  return(data)
}


#' Compute coverage data for a specific gene.
#'
#' This function computes data for a specific gene and remove intronic regions.
#'
#' @param coverageDataSet an objet of type \code{\link{coverageDataSet}}
#'        containing the coverage values for each sample.
#' @param name a character vector indicating the gene name.
#' @param introns an objet of type \code{\link{annotationsSet}}
#'        containing the annotations of the intronic regions.
#'        Note: The introns must have same name as the gene they
#'        are associated with.
#'
#' @return An objet of type geneCoverage containing the coverage values
#'         for the specified \code{gene} for all the samples in \code{coverageDataSet}.
#'
#' @import data.table
#' @export
#' @examples
#' data(rna_seq_signals)
#' data(introns)
#' gene_coverage(rna_seq_signals, "YER109C", introns)
gene_coverage <- function(coverageDataSet, name, introns = NULL) {
  if (!is.element(name, unique(coverageDataSet$data[, name]))) {
    stop(paste0("There is no data for ", name, "."))
  }

  if(!is.null(introns)){
    if(!is.element( "annotationsSet", class(introns))){
      introns <- as.annotationsSet(introns)
    }
  }

  gene_data <- coverageDataSet$data[name]

  gene_data$mean_type1 <- rowMeans(gene_data[, coverageDataSet$samples$type1, with = F])
  gene_data$mean_type2 <- rowMeans(gene_data[, coverageDataSet$samples$type2, with = F])

  ## remove introns regions
  if (!is.null(introns)) {
    clean_data <- .remove_introns(data = gene_data,
                                  name = name,
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
    gene_data <- lapply(gene_data[c(7:length(gene_data))], rev)
  }else{
    gene_data <- gene_data[c(7:length(gene_data))]
  }
  gene_data$original_pos <- clean_data$original_pos

  gene_data$gene_information <-
    as.data.frame(coverageDataSet$data[name][,1:6, with = F][1,])
  gene_data$introns <- clean_data$introns
  gene_data$samples <- coverageDataSet$samples

  class(gene_data) <- "geneCoverage"
  gene_data
}

.remove_introns <- function(data, name, introns) {

  original_pos <- data$pos

  if (is.element(name, introns$name)) {
    raw_intron <- introns[name]
    intron <- split(
      raw_intron,
      paste0("V", seq_len(nrow(raw_intron)))
    )
    gene_start <- data$start[1]

    intronic_regions = NULL
    intronic_regions <- lapply(intron, function(x) {
      if (x$end < gene_start | x$start > data$end[1]) {
        return(0)
      }
      (x$start - gene_start):(x$end - gene_start)
    })

    if (length(intronic_regions) == 1) {
      remove_rows <- base::unlist(intronic_regions)
      if (intronic_regions[[1]][1] != 0) {
        data <- data [!remove_rows]
        original_pos <- original_pos [-remove_rows]
      }
    }
  }else{
    data$pos <- seq_len(length(data$chrom))
    raw_intron <- NULL
  }
  list(original_pos = original_pos,
       data = data,
       introns = raw_intron)
}

## print generic function --------------------

#' @export
print.coverageDataSet <- function(x, ...) {
  l <- length(unique(x$data$name))
  if (l == 1) {
    cat ("\n\t\tcoverage values for ")
    cat (x$data$name[1])
    cat (" gene.\n\n")
  }else{
    cat ("\n\t\t\tcoverage dataset containing values for ")
    cat (l)
    cat (" genes.\n")
  }
  NextMethod("print", x$data)
}

#' @export
#' @importFrom utils head
print.geneCoverage <- function(x, ...) {
  cat ("\n\t\tCoverage data for gene: ")
  cat (x$gene_information$name)
  cat (". nrow = ")
  cat (length(x$pos))
  cat (". \n\n")

  l <- ifelse(is.null(x$introns), 2, 3)

  print(do.call(cbind, lapply(x[c(1:(length(x) - l))], head)))
}

# plot generic functions ----------------------

#' @S3method plot geneCoverage
#' @importFrom graphics abline axis legend lines par plot rect
#' @importFrom grDevices rgb
plot.geneCoverage <- function (x, cTSS = NULL, method = NULL, ...) {
  if (!is.null(cTSS)) {
    if (is.null(method)) {
      stop("please select a method you want to plot.\n
           (methodA, methodB, methodC, methodD, methodC_gaussian)")
    }
    method <- match.arg(method, c("methodC_gaussian", "methodA",
                                  "methodB", "methodC", "methodD"))
    cryptic_zones <- cTSS[method]
    }
  ylim <- max(x$mean_type1, x$mean_type2)
  # par(mar = c(5, 4, 1, 1),  las = 1)
  par(mar = c(3,2,0.5, 0.5), las = 1)
  plot(x = x$original_pos, y = x$mean_type1, ylim = c(0, ylim),
       xaxt = "n", bty = "n", ylab = "RNA-Seq signal (FPM)",
       xlab = "Distance from transcript start site (bp)", type = "l",
       cex.axis = 0.8, cex.lab = 0.8)
  l <- length(x$original_pos)
  axis(side = 1, at = c(seq(0, l, by = 500), l),
       labels = c(seq(0, l, by = 500), l), cex.axis = 0.8)
  lines(x = x$original_pos, y = x$mean_type2, lty = 2)
  samples <- c(strsplit(x$samples$type1[1], "_")[[1]][1],
               strsplit(x$samples$type2[1], "_")[[1]][1])

  legend("topleft", lty = c(1, 2), bty = "n", samples)
  if (!is.null(cTSS)) {
    ybottom = par()$usr[1]
    ytop = par()$usr[4]
    if (method == "methodC_gaussian") {
      invisible(mapply(function(mean, sd) {
        rect(xleft = mean - sd, xright = mean + sd, ybottom,
             ytop, col = grDevices::rgb(0, 0, 0, alpha = 0.3), border = F)
        abline(v = mean, col = "green")
      }, mean = cryptic_zones$methodC_gaussian$mean, sd = cryptic_zones$methodC_gaussian$sd))
    }
    else if (method == "methodD") {
      invisible(lapply(cryptic_zones, function(zone) {
        abline(v = zone[1]:zone[2], col = "darkgrey")
      }))
    }
    else {
      invisible(lapply(cryptic_zones, function(zone) {
        rect(xleft = zone[1], xright = zone[2], ybottom,
             ytop, col = rgb(0, 0, 0, alpha = 0.3), border = F)
      }))
    }
  }
}
