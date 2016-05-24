#' variantSetEnrichment
#'
#' This function will calculate matching random variant sets (MRVS) idential to AVS
#' @param avs A GRangesList object which is outputted by makeAVS function
#' @param mrvs A list of GRangesList objects outputted by makeMRVS function
#' @param regions A data frame. The data frame contains sample sheet identical to DiffBind or ChIPQC input sample sheets
#' @keywords VSE
#' @examples
#' \dontrun{
#' variantSetEnrichment(avs, mrvs, regions=samples)
#' #We have included the output from our example analysis with the package. \
#' You can load the example VariantSetEnrichment output by typing:
#' load(file.path(system.file("extdata", "vse_output.Rda", package="VSE")))
#' }
#' @importFrom car bcPower
#' @import GenomicRanges
#' @export
variantSetEnrichment <- function(avs, mrvs, regions){
  no_of_tags <- length(avs)
  no_of_beds <- length(regions$SampleID)
  intersect_matrix <- matrix(NA, nrow = no_of_beds, ncol = 1)
  row.names(intersect_matrix) <- regions$SampleID
  overlap <- matrix(NA, nrow = no_of_beds, ncol = no_of_tags)
  beds_list <- list()
  message("Loading regions")
  for (i in 1:no_of_beds){
    bed_path <- as.character(regions$Peaks[i])
    bed.gr <- bedToGRanges(bed_path)
    beds_list <- c(beds_list, bed.gr)
    overlap[i,] <- ifelse(countOverlaps(avs, bed.gr, ignore.strand=TRUE)>0,1,0)
    message(paste0(sum(overlap[i,]), " LD blocks intersect with ", regions$SampleID[i]))
  }
  intersect_matrix[,1] <- rowSums(overlap)
  null_intersects <- matrix(NA, nrow = no_of_beds, ncol = length(mrvs))
  row.names(null_intersects) <- regions$SampleID
  message("Tallying MRVS ", "\r", appendLF = FALSE)
  for (i in 1:length(mrvs)){
    flush.console()
    imod <- i %% length(mrvs)
    message(ifelse(imod %% 10 == 0, i, "."), "\r", appendLF = FALSE)
    overlap <- matrix(NA, nrow = no_of_beds, ncol = length(mrvs[[i]]))
    for (k in 1:length(beds_list)){
      overlap[k,] <- ifelse(countOverlaps(mrvs[[i]], beds_list[[k]], ignore.strand=TRUE)>0,1,0)
  #    print(overlap)
    }
    null_intersects[,i] <- rowSums(overlap)
#    print(null_intersects)
  };
  intersect_matrix <- cbind(intersect_matrix, null_intersects)

  #normalization
  normalityCutoff = 1
  matrix_norm <- t(apply(intersect_matrix, 1, function(x) x + runif(length(x),0,1)))
  ksvalues <- apply(matrix_norm, 1, function(x){
    kst <- ks.test(x, "pnorm", mean=mean(x), sd=sd(x), exact = TRUE)
    kst$p.value
  }    )
  message("\nNormalizing null distribution")
  for (i in 1:length(ksvalues)){
    if (ksvalues[i] < normalityCutoff){
      ks_lambda <- data.frame(lambda=seq(-2,2,0.1), ksp=rep(0, length(seq(-2,2,0.1))))
      for (j in 1:length(ks_lambda$lambda)){
        x <- car::bcPower(matrix_norm[i,], ks_lambda$lambda[j])
        kst <- ks.test(x, "pnorm", mean=mean(x), sd=sd(x), exact = TRUE)
        ks_lambda$ksp[j] <- kst$p.value
      }
      ks_lambda <- ks_lambda[order(ks_lambda$ksp, decreasing = TRUE),]
      matrix_norm[i,] <- car::bcPower(matrix_norm[i,], ks_lambda[1,1])
      ksvalues[i] <- ks_lambda[1,2]
    }
  }
  matrix_scaled <- t(apply(matrix_norm, 1, function(x) x <- (x - median(x))/ sd(x)))
  vse_matrices <- list(intersect_matrix, matrix_norm, matrix_scaled)
  return(vse_matrices)
}
