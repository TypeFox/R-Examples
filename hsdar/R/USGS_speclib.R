USGS_get_available_files <- function(url = NULL)
{
  if (!requireNamespace("RCurl", quietly = TRUE))
    stop("Library 'RCurl' is required to access USGS spectral library via ftp")
    
  if (is.null(url))
    url <- 'ftp://ftpext.cr.usgs.gov/pub/cr/co/denver/speclab/pub/spectral.library/splib06.library/ASCII/'
    
  if (substr(url, nchar(url), nchar(url)) != "/")
    url <- paste(url, "/", sep = "")
    
  filenames <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- strsplit(filenames, "\r*\n")[[1]]
  filenames.asc <- sapply(filenames, function(dir, url)
    {
      filenames <- RCurl::getURL(paste(url, dir, "/", sep = ""), ftp.use.epsv = FALSE, dirlistonly = TRUE)
      filenames <- strsplit(filenames, "\r*\n")[[1]]
      filenames.asc <- filenames[grep(".asc", filenames)]
      return(filenames.asc)
    }, url)
  attr(filenames.asc, "url") <- url
  return(filenames.asc)
}

USGS_retrieve_files <- function(avl = USGS_get_available_files(), pattern = NULL, retrieve = TRUE, loadAsSpeclib = TRUE, tol = 0.1)
{
  if (!requireNamespace("RCurl", quietly = TRUE))
    stop("Library 'RCurl' is required to access USGS spectral library via ftp")
  
  if (is.logical(loadAsSpeclib))
  {
    if (is.logical(retrieve))
    {
      if (!is.null(pattern))
      {
        for (i in 1:length(avl))
        {
          mat <- agrep(pattern, avl[[i]])
          avl[[i]] <- avl[[i]][mat]
        }
      }
    } else {
      avl <- retrieve
      retrieve <- TRUE    
    }
    if (!retrieve)
      return(avl)
    cat("Retrieve files...\n")
    files2read <- lapply(as.list(1:length(avl)), FUN = function(i, avl, url)
      {
        if (length(avl[[i]]) > 0)
        {
          return(sapply(avl[[i]], function(x, url, dir)
            {
              fi <- tempfile(fileext = ".asc")
              utils::download.file(paste(url, dir, "/", x, sep = ""), fi)
              return(fi)
            }, url, names(avl)[i]))
        } else {
          return(NA)
        }
      }, avl, attr(avl, "url"))
  } else {
    files2read <- loadAsSpeclib
    loadAsSpeclib <- TRUE
  }
  if (!loadAsSpeclib)
    return(files2read)
  cat(" done!\nRead files into speclib...")
  for (i in 1:length(avl))
  {
    if (length(avl[[i]]) > 0)
    {
      spec <- .read.USGS.asc(files2read[[i]], tol = tol)
      idSpeclib(spec) <- as.character(basename(avl[[i]]))
      if (!exists("ref"))
      {
        ref <- spec
      } else {
        ref <- .alignSpeclibs(ref, spec, tol = tol)
      }
    }
  }  
  notvalid <- apply(spectra(ref), 2, function(x) any(is.na(x)))
  ma <- as.vector(sapply(wavelength(ref)[notvalid], function(i) c(i-0.5,i+0.5)))
  if (length(ma) > 0)
    try(mask(ref) <- ma, silent = TRUE)
  cat(" done!\n")
  return(ref)
}

.read.USGS.asc <- function(filename, tol = 0.01)
{
  if (length(filename) > 1)
  {
    title <- as.character(1:length(filename))
    for (i in 1:length(filename))
    {
      if (i == 1)
      {
        ref <- .read.USGS.asc(filename[i])
        title[i] <- as.character(attribute(ref)$title[1])
      } else {
        dat <- .read.USGS.asc(filename[i])
        title[i] <- as.character(attribute(dat)$title[1])
        ref <- .alignSpeclibs(ref, dat, tol = tol)
      }      
    }
    attribute(ref) <- data.frame(title = title)
    idSpeclib(ref) <- as.character(basename(filename))
    return(ref)
  } else {
    fi <- file(filename, "r")
    titleline <- 0
    dataline <- 0
    start <- 0
    stop <- 0
    is.data <- FALSE
    nlines <- 0
    while (is.data != TRUE)
    {
      nlines <- nlines + 1
      line <- readLines(con = fi, n = 1)
      start <- start + 1
      if (substr(line,1,4) == "line")
      {
        tmp <- strsplit(line, " ")
        if (tmp[[1]][length(tmp[[1]])] == "title")
        {
          titleline <- as.numeric(tmp[[1]][2])
        } else {
          if (tmp[[1]][length(tmp[[1]])] == "history")
          {
            historyline <- tmp[[1]][2]
          } else {
            dataline <- as.numeric(tmp[[1]][2])            
          }
        }
      }
      if (nlines == titleline)
        title <- line
      if (nlines == dataline)
        is.data <- TRUE
    }    
    close(fi)
    dat <- read.table(filename, skip = dataline - 1, header = FALSE, dec = ".")
    dat[dat[,2] < 0, 2] <- NA
    dat <- speclib(spectra = dat[,2]*100, wavelength = dat[,1]*1000)
    attribute(dat) <- data.frame(title = title)
    return(dat)
  }
}

.alignSpeclibs <- function(x, y, tol = .01)
{
  al_attr <- TRUE
  test <- try(rbind(attribute(x), attribute(y)), silent = TRUE)
  if (inherits(test, "try-error"))
  {
    al_attr <- FALSE
    warning("Attribute information lost")
  }
  if (c("NONE", x@transformation)[length(x@transformation)+1] != c("NONE", y@transformation)[length(y@transformation)+1])
    stop("Transformation method between x and y differs")
  if (x@wlunit != y@wlunit)
    stop("Wavelength unit between x and y differs")
  if (any(c(length(x@rastermeta) > 0), length(y@rastermeta) > 0))
    warning("Rastermeta information will be lost")

  wl_1 <- round(wavelength(x)/tol, 0) * tol
  wl_2 <- round(wavelength(y)/tol, 0) * tol
  
  mat_1 <- match(wl_1, wl_2)
  mat_2 <- match(wl_2, wl_1)
  
  wl_3 <- c(wl_1[c(1:length(mat_1))[is.na(mat_1)]], wl_1[mat_2[!is.na(mat_2)]], wl_2[c(1:length(mat_2))[is.na(mat_2)]])
  
  mat_1 <- match(wl_1, wl_3, nomatch = 0)
  mat_2 <- match(wl_2, wl_3, nomatch = 0)
  
  spec <- matrix(NA, ncol = length(wl_3), nrow = nspectra(x) + nspectra(y))
  spec[1:nspectra(x), mat_1] <- spectra(x)
  spec[c((nspectra(x)+1):(nspectra(x) + nspectra(y))), mat_2] <- spectra(y)
  
  od <- order(wl_3)
  spec <- speclib(spec[,od], wavelength = wl_3[od])
  idSpeclib(spec) <- as.character(c(idSpeclib(x), idSpeclib(y)))
  
  if (al_attr)
    attribute(spec) <- rbind(attribute(x), attribute(y))
  
  return(spec)
}
  


