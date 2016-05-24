create_hapmap_reference <-
function(dir = getwd(),
                                    download_hapmap = FALSE,
                                    download_subset,
                                    hapmap_files = list.files(path = dir, pattern = "freqs_chr"),
                                    filename = "allele_reference_HapMap",
                                    save_txt = TRUE, save_rdata = !save_txt, return_reference = FALSE) {
  if(download_hapmap) {
    if(missing(download_subset)) stop("No subset selected. Options are: ASW, CEU, CHB, CHD, GIH, JPT, LWK, MEX, MKK, TSI, YRI")
    download_subset <- download_subset[1]
    if(!download_subset %in% c("ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI")) {
      stop("Subset argument unknown.  Options are: ASW, CEU, CHB, CHD, GIH, JPT, LWK, MEX, MKK, TSI, YRI") }
    
    print(paste("Downloading HapMap files - subset", download_subset), quote = FALSE)
    print(" - from: http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/2010-08_phaseII+III", quote = FALSE)
    print(" - Note: this data is downloaded from the HapMap website", quote = FALSE)
    print("   and subject to their terms and policies.", quote = FALSE)
    print("   See: http://hapmap.ncbi.nlm.nih.gov/datareleasepolicy.html", quote = FALSE)
    for(di in 1:25) {
      if(di < 23) { dn <- as.character(di)
      } else{
        if(di == 23L) dn <- "X"
        if(di == 24L) dn <- "Y"
        if(di == 25L) dn <- "M"
      }
      if(di != 1L) print(paste(" - downloading chromosome:", dn), quote = FALSE)
      flush.console()
      download.file(url = paste0("http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/2010-08_phaseII+III/",
                                 "allele_freqs_chr", dn, "_", download_subset, "_r28_nr.b36_fwd.txt.gz"),
                    destfile = paste0(dir, "/allele_freqs_chr", dn, "_", download_subset, "_r28_nr.b36_fwd.txt.gz"))
    }
    hapmap_files <- paste0("allele_freqs_chr", c(1:22, "X", "Y", "M"), "_", download_subset, "_r28_nr.b36_fwd.txt.gz")
    print("", quote = FALSE)
    print("Processing files", quote = FALSE)
  } else {
    print("Reading files:", quote = FALSE)
    print(hapmap_files, quote = FALSE)   
  }

  flush.console()
  fileN <- length(hapmap_files)
  if(fileN == 0L) stop("No hapmap_files specified")
  for(fi in 1:fileN) {		
    if(!file.exists(paste(dir, hapmap_files[fi], sep = "/"))) stop(paste("cannot find", hapmap_files[fi]))
    chrmap <- read.table(paste(dir, hapmap_files[fi], sep = "/"), comment.char = "", header = TRUE, stringsAsFactors = FALSE)[ , c(1:3, 14, 11, 15, 12)]
    
    chr <- substr(chrmap$chrom[1], 4L, nchar(chrmap$chrom[1]) )
    if(chr == "X") chr <- 23L
    if(chr == "Y") chr <- 24L
    if(chr == "XY")chr <- 25L
    if(chr == "M") chr <- 26L
    chrmap$chrom <- as.integer(chr)
    
    flip <- which(chrmap$otherallele_freq > 0.5 & chrmap$otherallele_freq > chrmap$refallele_freq)
    if(length(flip) > 0L) {
      tempallele				<- chrmap$refallele[flip]
      chrmap$refallele[flip]		<- chrmap$otherallele[flip]
      chrmap$otherallele[flip]	<- tempallele
      tempfrq				<- chrmap$refallele_freq[flip]
      chrmap$refallele_freq[flip]	<- chrmap$otherallele_freq[flip]
      chrmap$otherallele_freq[flip] <- tempfrq
    }
    
    print(paste0("Finished chr ", chr, " (", fi, " out of ", fileN, ")  - SNPs: ", nrow(chrmap), " - flipped: ", length(flip)), quote = FALSE)
    flush.console()
    allele_ref_std <- if(fi == 1L) chrmap else rbind(allele_ref_std, chrmap)
  }
  colnames(allele_ref_std) <- c("SNP", "CHR", "POS", "MINOR", "MAJOR", "MAF", "other_AF")
  
  allele_ref_std <- allele_ref_std[order(allele_ref_std$CHR, na.last = FALSE), ]
  
  if(sum(duplicated(allele_ref_std$SNP)) > 0L) {
    dupli_names <- unique(allele_ref_std$SNP[duplicated(allele_ref_std$SNP)])
    duplis <- which(allele_ref_std$SNP %in% dupli_names)
    print("WARNING: duplicate SNPids encountered in reference dataset - SNPs removed", quote = FALSE)
    write.table(allele_ref_std[duplis, ], paste0(dir, "/", filename, "_removed_duplicates.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  } else { duplis <- NULL }
  
  monomorph <- which(!allele_ref_std$MINOR %in% c("A", "T", "C", "G") | !allele_ref_std$MAJOR %in% c("A", "T", "C", "G"))
  if(length(monomorph) > 0L) {
    print("WARNING: invalid / missing alleles encountered in reference dataset - SNPs removed", quote = FALSE)
    write.table(allele_ref_std[monomorph,], paste0(dir, "/", filename, "_removed_bad_alleles.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  }
  
  bad_freq <- which(allele_ref_std$MAF + allele_ref_std$other_AF != 1)
  if(length(bad_freq ) > 0L) {
    print("WARNING: allele frequencies do not add up - SNPs removed", quote = FALSE)
    write.table(allele_ref_std[bad_freq,], paste0(dir, "/", filename, "_bad_frequency.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  }
  
  remove_list <- unique(c(duplis, monomorph, bad_freq))
  if(length(remove_list) > 0L) {
    if(nrow(allele_ref_std) == length(remove_list)) stop("No SNPs remaining in dataset - the output files should indicate why all SNPs are excluded")
    allele_ref_std <- allele_ref_std[-remove_list, 1:6]
  }
  
  allele_ref_std$SOURCE <- as.factor(filename)
  allele_ref_std$DATE_ADDED <- as.factor(date())

  if(download_hapmap) {
   print("", quote = FALSE)
   print("The downloaded HapMap files have been saved in:", quote = FALSE)
   print(paste("  ", dir), quote = FALSE)
   print("   They are no longer necessary and can be deleted.", quote = FALSE)
   print(" - Note: these files are downloaded from the HapMap website", quote = FALSE)
   print("   and subject to their terms and policies.", quote = FALSE)
   print("   See: http://hapmap.ncbi.nlm.nih.gov/datareleasepolicy.html", quote = FALSE)
  }
  
  if(save_txt) write.table(allele_ref_std, paste0(dir, "/", filename, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  if(save_rdata | return_reference) {
    allele_ref_std$MINOR <- as.factor(allele_ref_std$MINOR)
    allele_ref_std$MAJOR <- as.factor(allele_ref_std$MAJOR)
    if(save_rdata) save(allele_ref_std, file = paste0(dir, "/", filename, ".RData"))
    if(return_reference) return(allele_ref_std)
  }
  
  return(invisible())
}
