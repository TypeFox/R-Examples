################################################################################
# haploDiv: calculate Weir & Cockerham's Fst from haploid genotypes            #
################################################################################
# haploDiv function for calculating various statistics from haploid data
# try diploidization first
#' @export
haploDiv <- function(infile = NULL, outfile = NULL, pairwise = FALSE, 
                     boots = 0){
  if(boots != 0){
    bs_pairwise <- TRUE
    para <- TRUE
  } else {
    bs_pairwise <- FALSE
    para <- FALSE
  }
  haploFileReader <- function(x){
    fileReader <- function(infile){
      if (typeof(infile) == "list") {
        return(infile)
      } else if (typeof(infile) == "character") {
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if (ext == "arp") {
          convRes <- arp2gen(infile)
          if (!is.null(convRes)) {
            cat("Arlequin file converted to genepop format! \n")
            infile <- paste(flForm[1], ".gen", sep = "")
          } else {
            infile <- paste(flForm[1], ".gen", sep = "")
          }
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
          locs <- strsplit(dat[2], split = "\\s+")[[1]]
          if(length(locs != 1)){
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
        
        
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if (popLoc[1] == 3) {
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:length(dat)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x) {
          x <- unlist(strsplit(x, split = "\\s+"))
          if (is.element("", x)) {
            x <- x[-(which(x == ""))]
          }
          if (is.element(",", x)) {
            x <- x[-(which(x == ","))]
          }
          if (length(x) != 1 && length(x) != no_col) {
            x <- paste(x, collapse = "")
          }
          if (length(x) < no_col) {
            tabs <- paste(rep(NA, (no_col - length(x))), 
                          sep = "\t", collapse = "\t")
            line <- paste(x, tabs, sep = "\t")
            line <- unlist(strsplit(line, split = "\t"))
            return(line)
          } else {
            return(x)
          }
        })
      }
      out <- as.data.frame(t(dat1))
      rownames(out) <- NULL
      return(out)
    }
    z <- as.matrix(fileReader(x))
    nloc <- ncol(z)-1
    popStrt <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(z[,1])) + 1
    popEnd <- c(popStrt[-1] - 2, nrow(z))
    gp <- nchar(as.character(z[popStrt[1], 2]))
    diploFun <- function(strt, end, dat, gp){
      tst <- t(sapply(strt:end, function(i){
        if(gp == 2){
          paste(sprintf("%02g", as.numeric(dat[i, -1])), 
                sprintf("%02g", as.numeric(dat[i, -1])), sep = "")
        } else if(gp == 3){
          paste(sprintf("%03g", as.numeric(dat[i, -1])), 
                sprintf("%03g", as.numeric(dat[i, -1])), sep = "")
        } else {
          cat("There is a problem with your input file!")
        }
      }))
      return(tst)
    }
    
    diploidGeno <- mapply(diploFun, strt = popStrt, end = popEnd, 
                          MoreArgs = list(dat = z, gp = gp), SIMPLIFY = FALSE)
    dat <- as.matrix(z)
    for(i in 1:length(popEnd)){
      dat[popStrt[i]:popEnd[i], -1] <- diploidGeno[[i]]
    }
    list(data = as.data.frame(dat),
         gp = gp)
  }
  # read the file and diploidize
  dat <- haploFileReader(infile)
  out <- diveRsity::fastDivPart(infile = dat$data, outfile = NULL, 
                                gp = dat$gp, pairwise = pairwise, 
                                fst = TRUE, boots = boots,
                                bs_pairwise = bs_pairwise, para = para)
  if(pairwise && boots > 0){
    output <- list(locus = out$estimate[-nrow(out$estimate), "Fst_WC"],
                   overall = out$estimate[nrow(out$estimate), "Fst_WC"],
                   pairwise = out$pairwise$thetaWC,
                   bs_pairwise = out$bs_pairwise$thetaWC) 
  } else if(pairwise && boots == 0L){
    output <- list(locus = out$estimate[-nrow(out$estimate), "Fst_WC"],
                   overall = out$estimate[nrow(out$estimate), "Fst_WC"],
                   pairwise = out$pairwise$thetaWC) 
  } else{
    output <- list(locus = out$estimate[-nrow(out$estimate), "Fst_WC"],
                   overall = out$estimate[nrow(out$estimate), "Fst_WC"])
  }
  
  # write reuslts to file
  if(!is.null(outfile) && pairwise && bs_pairwise){
    pwmat <- round(output$pairwise, 4)
    pwmat[is.na(pwmat)] <- ""
    idx <- 1:ncol(pwmat)
    for(i in 1:length(idx)){
      pwmat[idx[i], idx[i]] <- "--"
    }
    # pairwise matrix
    fl <- file(paste(outfile, "-[pwMatrix].txt", sep = ""), "w")
    cat("Pairwise Fst (Weir & Cockerham, 1984)", "\n", sep = "\t", file = fl)
    cat("", "\n", sep = "\t", file = fl)
    cat(c("", colnames(output$pairwise)), "\n", sep = "\t", file = fl)
    for(i in 1:nrow(output$pairwise)){
      cat(c(colnames(output$pairwise)[i], pwmat[i, ]), 
          "\n", sep = "\t", file = fl)
    }
    close(fl)
    # Pairwise cis
    fl <- file(paste(outfile, "-[pwBootstrap].txt", sep = ""), "w")
    cat("Bootstrapped 95% Confidence intervals for Weir & Cockerham's (1984) Fst",
        "\n", sep = "\t", file = fl)
    cat("", "\n", file = fl)
    cat(c("","actual", "mean", "BC_mean", "lower", "upper", 
          "BC_lower", "BC_upper"), "\n", sep = "\t", file = fl)
    for(i in 1:nrow(output$bs_pairwise)){
      cat(c(rownames(output$bs_pairwise)[i], round(output$bs_pairwise[i, ], 4)), 
          "\n", sep = "\t", file = fl)
    }
    close(fl)
  } else if(!is.null(outfile) && pairwise && !bs_pairwise){
    # pairwise matrix
    pwmat <- round(output$pairwise, 4)
    pwmat[is.na(pwmat)] <- ""
    idx <- 1:ncol(pwmat)
    for(i in 1:length(idx)){
      pwmat[idx[i], idx[i]] <- "--"
    }
    fl <- file(paste(outfile, "-[pwMatrix].txt", sep = ""), "w")
    cat("Pairwise Fst (Weir & Cockerham, 1984)", "\n", sep = "\t", file = fl)
    cat("", "\n", sep = "\t", file = fl)
    cat(c("", colnames(output$pairwise)), "\n", sep = "\t", file = fl)
    for(i in 1:nrow(output$pairwise)){
      cat(c(colnames(output$pairwise)[i], pwmat[i, ]), 
          "\n", sep = "\t", file = fl)
    }
    close(fl)
  }
  return(output)
}
################################################################################
# haploDiv end                                                                 #
################################################################################