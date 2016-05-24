################################################################################
# bigDivPart - a wrapper function for the calculation of diff stats
################################################################################
#' @export
bigDivPart <- function(infile = NULL, outfile = NULL, WC_Fst = FALSE,
                       format = NULL){
  .Deprecated(new = "diffCalc", msg = "This function is no longer in use. Please use 'diffCalc' instead, \nSee ?diffCalc for usage details.", 
              old = "bigDivPart")
}

# old function definition
# bigDivPart <- function(infile = NULL, outfile = NULL, WC_Fst = FALSE,
#                        format = NULL){
#   fstat = WC_Fst
#   on = outfile
#   if (!is.null(on) && format != "txt" && format != "xlsx") {
#     stop("Please provide a valid output file format")
#   }
#   fastScan <- function(fname) {
#     s <- file.info(fname)$size
#     buf <- readChar(fname, s, useBytes = TRUE)
#     # replace Mac encoded line endings
#     if(length(grep("\r", buf)) != 0L){
#       buf <- gsub("\r", "\n", buf)
#       buf <- gsub("\n\n", "\n", buf)
#     }
#     return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
#   }
#   dat <- fastScan(fname = infile)
#   if (length(strsplit(dat[length(dat)], split = "\\s+")[[1]]) == 
#         1) {
#     dat <- dat[-(length(dat))]
#   }
#   rm(fastScan)
#   z <- gc()
#   rm(z)
#   popLocation <- grep("^([[:space:]]*)POP([[:space:]]*)$", 
#                       toupper(dat))
#   pop_pos <- c(popLocation, (length(dat) + 1))
#   loci_names <- as.vector(sapply(dat[2:(pop_pos[1] - 1)], function(x) {
#     gsub(pattern = "\\s+", replacement = "", x)
#   }))
#   popSizes <- NULL
#   for (i in 1:(length(pop_pos) - 1)) {
#     popSizes[i] <- length((pop_pos[i] + 1):(pop_pos[(i + 
#                                                        1)] - 1))
#   }
#   pops <- dat[-(c(1:(popLocation[1] - 1), popLocation))]
#   popList <- lapply(seq_along(popSizes), function(i) {
#     if (i == 1) {
#       indx <- 1:popSizes[i]
#     } else {
#       indx <- (sum(popSizes[1:(i - 1)]) + 1):((sum(popSizes[1:(i - 1)])) +
#                                                 popSizes[i])
#     }
#     return(pops[indx])
#   })
#   npops <- length(popList)
#   nloci <- length(loci_names)
#   pop_sizes <- popSizes
#   rm(dat, pops)
#   z <- gc(reset = TRUE)
#   rm(z)
#   testStr <- strsplit(popList[[1]][1], split = "\\s+")[[1]]
#   gpEst <- sapply(testStr, function(x) {
#     if (is.character(x)) {
#       nchar(x)/2
#     } else {
#       NA
#     }
#   })
#   rm(testStr)
#   gp <- as.numeric(names(sort(-table(gpEst)))[1])
#   prePopList <- lapply(popList, function(x) {
#     y <- array(data = NA, dim = c(length(x), (nloci + 1), 2))
#     colnames(y) <- c("ind", loci_names)
#     for (j in 1:length(x)) {
#       data <- strsplit(x[j], split = "\\s+")[[1]]
#       if (data[2] == ",") {
#         data <- data[-2]
#       }
#       data[data == "NANA"] <- NA
#       data[data == "0"] <- NA
#       data[data == "000000"] <- NA
#       data[data == "999999"] <- NA
#       data[data == "-9-9"] <- NA
#       data[data == "0000"] <- NA
#       y[j, 2:(nloci + 1), 1] <- substr(data[2:(nloci + 1)], 1, gp)
#       y[j, 2:(nloci + 1), 2] <- substr(data[2:(nloci + 1)], gp + 1, gp * 2)
#       y[j, 1, 1] <- data[1]
#       y[j, 1, 2] <- data[1]
#     }
#     return(y)
#   })
#   rm(popList)
#   ind_names <- lapply(prePopList, function(x) {
#     return(x[, 1, 1])
#   })
#   pop_names <- sapply(ind_names, function(x) {
#     return(x[1])
#   })
#   nb <- bigPreDiv(prePopList, FALSE, nloci, npops, popSizes, 
#                   fstat)
#   stdOut <- data.frame(loci = c(loci_names, "Global"), 
#                        H_st = c(nb$hst, NA), 
#                        D_st = c(nb$dst, NA), 
#                        G_st = c(nb$gst, nb$gst_all), 
#                        G_hed_st = c(nb$gst_hedrick, nb$gst_all_hedrick), 
#                        D_Jost = c(nb$djost, nb$djost_all))
#   if (fstat) {
#     estOut <- data.frame(loci = c(loci_names, "Global"), 
#                          Harmonic_N = c(nb$locus_harmonic_N, NA), 
#                          H_st_est = c(nb$hst_est, NA), 
#                          D_st_est = c(nb$dst_est, NA), 
#                          G_st_est = c(nb$gst_est, nb$gst_est_all), 
#                          G_hed_st = c(nb$gst_est_hedrick, 
#                                       nb$gst_est_all_hedrick), 
#                          D_Jost = c(nb$djost_est, nb$djost_est_all), 
#                          Fst_WC = nb$fstats[, 1], Fit_WC = nb$fstats[, 2])
#   } else {
#     estOut <- data.frame(loci = c(loci_names, "Global"), 
#                          Harmonic_N = c(nb$locus_harmonic_N, NA), 
#                          H_st_est = c(nb$hst_est, NA), 
#                          D_st_est = c(nb$dst_est, NA), 
#                          G_st_est = c(nb$gst_est, nb$gst_est_all), 
#                          G_hed_st = c(nb$gst_est_hedrick, 
#                                       nb$gst_est_all_hedrick), 
#                          D_Jost = c(nb$djost_est, nb$djost_est_all))
#   }
#   if (!is.null(on)) {
#     suppressWarnings(dir.create(path = paste(getwd(), "/", 
#                                              on, "-[diveRsity]", "/", 
#                                              sep = "")))
#     of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
#   }
#   write_res <- is.element("xlsx", installed.packages()[, 1])
#   if (!is.null(on)) {
#     if (write_res && format == "xlsx") {
#       require("xlsx")
#       write.xlsx(stdOut, file = paste(of, "[bigDivPart].xlsx", 
#                                       sep = ""), 
#                  sheetName = "Standard_stats", col.names = TRUE, 
#                  row.names = FALSE, append = FALSE)
#       
#       write.xlsx(estOut, file = paste(of, "[bigDivPart].xlsx", 
#                                       sep = ""), 
#                  sheetName = "Estimated_stats", col.names = TRUE, 
#                  row.names = FALSE, append = TRUE)
#     } else {
#       std <- file(paste(of, "Standard-stats[bigDivPart].txt", 
#                         sep = ""), "w")
#       cat(paste(colnames(stdOut), sep = ""), "\n", sep = "\t", 
#           file = std)
#       stdOut <- as.matrix(stdOut)
#       for (i in 1:nrow(stdOut)) {
#         cat(stdOut[i, ], "\n", file = std, sep = "\t")
#       }
#       close(std)
#       est <- file(paste(of, "Estimated-stats[bigDivPart].txt", 
#                         sep = ""), "w")
#       cat(paste(colnames(estOut), sep = ""), "\n", sep = "\t", 
#           file = est)
#       estOut <- as.matrix(estOut)
#       for (i in 1:nrow(estOut)) {
#         cat(estOut[i, ], "\n", file = est, sep = "\t")
#       }
#       close(est)
#     }
#   }
#   list(standard = stdOut, estimates = estOut)
#}
################################################################################
#  END - bigDivPart
################################################################################