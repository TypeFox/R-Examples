atcg1234 <- function(data, ploidy=2, format="ATCG", maf=0, multi=TRUE){
  #### apply with progress bar ######
  apply_pb <- function(X, MARGIN, FUN, ...)
  {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTxtProgressBar(get("pb", envir= env),
                        curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  # tells you which markers have double letter code, i.e. TT instead of T
  # 1: has only one letter
  # 0: has two letters
  user.code <- apply(data[,c(1:(round(dim(data)[2]/20)))], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
  AA <- sum(user.code, na.rm = TRUE)/length(user.code)
  if(AA > .9){
    rrn <- rownames(data)
    ##### START GBS.TO.BISNP DATA ######
    gbs.to.bisnp <- function(x) {
      y <- rep(NA,length(x))
      y[which(x=="A")] <- "AA"
      y[which(x=="T")] <- "TT"
      y[which(x=="C")] <- "CC"
      y[which(x=="G")] <- "GG"
      y[which(x=="R")] <- "AG"
      y[which(x=="Y")] <- "CT"
      y[which(x=="S")] <- "CG"
      y[which(x=="W")] <- "AT"
      y[which(x=="K")] <- "GT"
      y[which(x=="M")] <- "AC"
      y[which(x=="+")] <- "NN"
      y[which(x=="0")] <- "NN"
      y[which(x=="-")] <- NA
      y[which(x=="N")] <- NA
      return(y)
    }
    ##### END GBS.TO.BISNP DATA ######
    cat("Converting GBS single-letter to biallelic code\n")
    data <- apply_pb(data, 2,gbs.to.bisnp)
    rownames(data) <- rrn
    data <- as.data.frame(data)
  }
  #### apply with progress bar ######
  s1 <- rownames(data)
  s2 <- colnames(data)
  data <- as.data.frame(t(data))
  rownames(data) <- s2
  colnames(data) <- s1
  bases <- c("A", "C", "G", "T","l","m","n","p","h","k","-","+")
  ## get reference allele function
  get.ref <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ATCG") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        ref.alt <- (bases[which(ans == 1)])[1:2]
        #stop("Error in genotype matrix: More than 2 alleles")
      }
      if (sum(ans) == 2) {
        ref.alt <- bases[which(ans == 1)]
      }
      if (sum(ans) == 1) {
        ref.alt <- c(bases[which(ans == 1)], NA)
      }
    }
    return(ref.alt)
  }
  
  get.multi <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ATCG") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        ref.alt <- TRUE
      }
      if (sum(ans) == 2) {
        ref.alt <- FALSE
      }
      if (sum(ans) == 1) {
        ref.alt <- FALSE
      }
    }
    return(ref.alt)
  }
  
  ####################################
  ## convert to matrix format
  ####################################
  markers <- as.matrix(data)
  ####################################
  # get reference alleles
  ####################################
  cat("Obtaining reference alleles\n")
  tmp <- apply_pb(markers, 1, get.ref, format=format)
 
  if(multi){ # if markers with multiple alleles should be removed
    cat("Checking for markers with more than 2 alleles. If found will be removed.\n")
    tmpo <- apply_pb(markers, 1, get.multi, format = format)
    multi.allelic <- which(!tmpo) # good markers
    markers <- markers[multi.allelic,]
    tmp <- tmp[, multi.allelic]
  }
  
  Ref <- tmp[1, ]
  Alt <- tmp[2, ]
  ####################################
  ## bind reference allele and markers and convert to numeric format based on the 
  # reference/alternate allele found
  ####################################
  cat("Converting to numeric format\n")
  M <- apply_pb(cbind(Ref, markers), 1, function(x) {
    y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
    ans <- as.integer(lapply(y, function(z) {
      ifelse(z[1] < 0, ploidy, ploidy - length(z))
    }))
    return(ans)
  })
  gid.geno <- s1 #colnames(geno)
  rownames(M) <- gid.geno
  ####################################
  # identify bad markers
  ####################################
  bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
  if (bad > 0) {
    stop("Invalid marker calls.")
  }
  ####################################
  # by column or markers calculate MAF
  ####################################
  cat("Calculating minor allele frequency (MAF)\n")
  MAF <- apply_pb(M, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  ####################################
  # which markers have MAF > 0, JUST GET THOSE
  ####################################
  polymorphic <- which(MAF > maf)
  M <- M[, polymorphic]
  ####################################
  # function to impute markers with the mode
  ####################################
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  # time to impute
  missing <- which(is.na(M))
  if (length(missing) > 0) {
    cat("Imputing missing data with mode \n")
    M <- apply_pb(M, 2, impute.mode)
  }
  if(ploidy == 2){
    M <- M - 1
  }
  #rownames(M) <- rownames(data)
  ####################################
  return(M)
}