################################################################################
# rpg: a new faster, memory efficient function for reading genepop files       #
################################################################################
#' New readgenepop format
#' 
#' Kevin Keenan 2014
rgp <- function(infile){
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    # replace Mac encoded line endings
    if(length(grep("\r", buf)) != 0L){
      buf <- gsub("\r", "\n", buf)
      buf <- gsub("\n\n", "\n", buf)
    }
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  if(is.list(infile)){
    infile <- as.matrix(infile)
    dat <- apply(infile, 1, function(x){
      x <- x[!is.na(x)]
      return(paste(x, collapse = "\t"))
    })
    #dat <- c(paste(colnames(infile), collapse = "\t"), dat)
  } else {
    dat <- fastScan(infile)
    # strip whitespace from the beginning an end of lines
    dat <- sapply(dat, function(x){
      sub("^\\s+", "", x)
    })
    dat <- sapply(dat, function(x){
      return(sub("\\s+$", "", x))
    })
    names(dat) <- NULL
  }
  popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
  if(popLoc[1] == 3){
    if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
      locs <- strsplit(dat[2], split = "\\s+")[[1]]
      if(length(locs) == 1){
        locs <- strsplit(dat[2], split = ",")[[1]]
      }
      locs <- as.character(sapply(locs, function(x){
        x <- strsplit(x, split = "")[[1]]
        if(is.element(",", x)){
          x <- x[-(which(x == ","))]
        }
        return(paste(x, collapse = ""))
      }))
      dat <- c(dat[1], locs, dat[-(1:2)])
    }
  } else {
    locs <- as.character(dat[2:(popLoc[1]-1)])
  }
  # strip whitespace from locus names
  locs <- as.character(sapply(locs, function(x){
    return(strsplit(x, split = "\\s+")[[1]][1])
  }))
  # npops
  popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
  npops <- length(popLoc)
  no_col <- length(locs)+1
  nloci <- length(locs)
  # get genotypes
  strt <- popLoc + 1
  ends <- c(popLoc[-1] - 1, length(dat))
  genoRet <- function(strt, ends, x){
    out <- strsplit(x[strt:ends], split = "\\s+")
    x <- do.call("rbind", c(out, deparse.level = 0))
    if(round(mean(nchar(x[,2]))) == 1L){
      x[,1] <- paste(x[,1], x[,2], sep = "")
      x <- x[,(-2)]
    }
    x[x == "-9"] <- NA
    x[x == "0000"] <- NA
    x[x == "000000"] <- NA
    # output
    list(ls = x[,(-1)],
         nms = as.vector(x[,1]))
  }
  genos <- mapply(genoRet, strt = strt, ends = ends, 
                  MoreArgs = list(x = dat), SIMPLIFY = FALSE)
  indNames <- lapply(genos, "[[", 2)
  #indNames <- do.call("c", indNames)
  genos <- lapply(genos, "[[", 1)
  # detect genepop format
  # check for loci with all missing data before calculating gp
  badLoc <- apply(genos[[1]], 2, function(x){
    sum(is.na(x)) == length(x)
  })
  badLoc <- which(!badLoc)
  gp <- round(mean(nchar(na.omit(genos[[badLoc[1]]][,1]))/2))
  # convert genotypes to arrays
  genos <- lapply(genos, function(x){
    al1 <- substr(x, 1, gp)
    al2 <- substr(x, (gp+1), (gp*2))
    out <- array(NA, dim = c(nrow(x), ncol(x), 2))
    out[,,1] <- al1
    out[,,2] <- al2
    return(out)
  })
  
  # calculate allele frequencies, obs alleles, popSizes
  # define function
  statFun <- function(x, cl = NULL){
    # if(!is.null(cl)){
    # tab <- parLapply(cl, 1:dim(x)[2], function(i){return(table(x[,i,]))})
    # } else {
    #tab <- lapply(1:dim(x)[2], function(i){return(table(x[,i,]))})
    #}
    popSizes <- apply(x, 2, function(y){
      length(na.omit(y[,1])) * 2
    })
    af <- lapply(1:dim(x)[2], function(i){
      y <- as.vector(na.omit(x[,i,]))
      nms <- unique(y)[order(unique(y))]
      ot <- myTab(y)
      names(ot) <- nms
      return(ot)
    })
    popSizes <- popSizes/2
    list(af = af, ps = popSizes)
  }
  # rearrange data by loci
  check <- function(args, gp){
    #args <- list(...)
    npops <- length(args)
    pad <- paste("%0", gp, "g", sep = "")
    rnames <- sprintf(pad, 
                      unique(sort(as.numeric(unlist(lapply(args,
                                                           names))))))
    
    out <- matrix(0, nrow = length(rnames), ncol = npops)
    rownames(out) <- as.character(rnames)
    for(i in 1:npops){
      out[match(names(args[[i]]), rownames(out)),i] <- as.numeric(args[[i]])
    }
    return(out)
  }
  # calculate stats
  obsAllSize <- lapply(genos, statFun)
  # get individual stats
  af <- lapply(obsAllSize, function(x){
    out <- x$af
    x$af <- NULL
    return(out)
  })
  #obs <- lapply(obsAllSize, function(x){
  #  return(x$obs)
  #})
  ps <- lapply(obsAllSize, function(x){
    out <- x$ps
    x$ps <- NULL
    return(out)
  })
  af <- lapply(1:(nloci), function(i){
    return(lapply(af, "[[", i))
  })
  #obs <- lapply(1:(nloci), function(i){
  #  return(lapply(obs, "[[", i))
  #})
  ps <- lapply(1:(nloci), function(i){
    return(sapply(ps, "[", i))
  })
  af <- lapply(af, check, gp = gp)
  # names(af) <- locs
  #obs <- lapply(obs, check)
  gc()
  list(af = af, genos = genos, ps = ps, gp = gp,
       indnms = indNames, locs = locs)
}
################################################################################
# end rpg                                                                      #
################################################################################