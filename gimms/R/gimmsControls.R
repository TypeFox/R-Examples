### download function ----------------------------------------------------------

downloader <- function(x, dsn = getwd(), overwrite = FALSE, quiet = TRUE,
                       mode = "wb", cores = 1L, ...) {

  ### single core

  if (cores == 1L) {

    ## download files one after another
    for (i in x) {
      destfile <- paste0(dsn, "/", basename(i))
      if (file.exists(destfile) & !overwrite) {
        if (!quiet)
          cat("File", destfile, "already exists in destination folder. Proceeding to next file ...\n")
      } else {
        try(download.file(i, destfile = destfile, mode = mode,
                          quiet = quiet, ...), silent = TRUE)
      }
    }


    ### multi-core

  } else {

    ## initialize cluster
    cl <- parallel::makePSOCKcluster(cores)

    ## export required variables
    parallel::clusterExport(cl, c("x", "cores", "dsn", "overwrite", "quiet",
                                  "mode"), envir = environment())

    ## download files in parallel
    parallel::parLapply(cl, x, function(i) {
      destfile <- paste0(dsn, "/", basename(i))
      if (file.exists(destfile) & !overwrite) {
        if (!quiet)
          cat("File", destfile, "already exists in destination folder. Proceeding to next file ...\n")
      } else {
        try(download.file(i, destfile = destfile, mode = mode,
                          quiet = quiet, ...), silent = TRUE)
      }
    })

    ## deregister parallel backend
    parallel::stopCluster(cl)
  }

  ## return downloaded files
  gimms_out <- paste(dsn, basename(x), sep = "/")
  return(gimms_out)
}


### create gimms-specific envi header file -------------------------------------

createHeader <- function(file, header) {

  ## location of header file
  file_hdr <- paste0(file, ".hdr")

  ## default content of gimms ndvi3g-related header file
  if (missing(header))
    header <- paste("ENVI",
                    "description = { R-language data }",
                    "samples = 2160",
                    "lines = 4320",
                    "bands = 1",
                    "data type = 2",
                    "header offset = 0",
                    "interleave = bsq",
                    "sensor type = AVHRR",
                    "byte order = 1", sep = "\n")

  ## write and return file
  writeLines(header, file_hdr)
  return(file_hdr)
}


### check desired number of cores ----------------------------------------------

checkCores <- function(cores) {

  ## available cores
  cores_avl <- parallel::detectCores()

  ## resize if 'cores' exceeds number of available cores
  if (cores > cores_avl) {
    cores <- cores_avl - 1
    warning("Desired number of cores is invalid. Resizing parallel cluster to ", cores, " cores.")
  }

  return(cores)
}
