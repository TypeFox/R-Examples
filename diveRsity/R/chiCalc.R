#' chiCalc function
#' Fisher's exact testing of sample independence from genotype data.
#' Allows the calculation of pairwise differences
#' 
#' Kevin Keenan, QUB, 2014
#' 
#' @export
chiCalc <- function(infile = NULL, outfile = NULL, pairwise = FALSE,
                     mcRep = 2000){ 
  dat <- rgp(infile)
  #pairwise = T
  #mcRep = 2000
  #outfile <- "test"
  # calulate observed genotypes
  obs <- lapply(dat$genos, function(x){
    lapply(1:dim(x)[2], function(i){
      y <- x[,i,]
      ip <- apply(y, 1, paste, collapse = "")
      ip[ip == "NANA"] <- NA
      ip <- as.vector(na.omit(ip))
      Tab(ip)
    })
  })
  obs <- lapply(1:length(obs[[1]]), function(i) lapply(obs, "[[", i))
  if(pairwise){
    pw <- combn(length(dat$genos), 2)
  }
  # create a chi-dif function
  chiDif <- function(..., mcRep){
    ip <- (...)
    nms <- unique(names(unlist(ip)))
    if(length(nms) == 1L){
      nms <- c(nms, "NA")
    }
    mat <- matrix(0, ncol = length(nms), nrow = length(ip))
    colnames(mat) <- nms
    for(i in 1:length(ip)){
      mat[i, names(ip[[i]])] <- ip[[i]]
    }
    fisher.test(mat, simulate.p.value = TRUE, B = mcRep)$p.value
  }
  # calculate overall chisq stats
  glb <- sapply(obs, chiDif, mcRep = mcRep)
  chiAll <- function(x){
    df <- 2*length(x)
    chi <- -2*sum(log(x), na.rm = TRUE)
    p <- pchisq(chi, df = df, lower.tail = FALSE)
    list(p = p, chisq = chi)
  }
  all <- chiAll(glb)$p
  glb <- data.frame(loci = c(dat$locs, "Overall"), 
                    p.value = round(c(glb, all), 4))
  if(pairwise){
    # calculate pairwise differentiation
    pwDif <- lapply(1:ncol(pw), function(i){
      pwobs <- lapply(obs, "[", pw[,i])
      pwp <- sapply(pwobs, chiDif, mcRep = mcRep)
      all <- chiAll(pwp)$p
      return(data.frame(loci = c(dat$locs, "Overall"),
                        p.value = round(c(pwp, all), 4)))
    })
    # pairwise names
    nms <- sapply(dat$indnms, "[", 1)
    pwnms <- paste(nms[pw[1,]], nms[pw[2,]], sep = " vs ")
    ovrallpw <- sapply(pwDif, function(x){
      return(x$p.value[x$loci == "Overall"])
    })
    pwMat <- matrix(NA, ncol = length(nms), nrow = length(nms),
                    dimnames = list(nms, nms))
    for(i in 1:ncol(pw)){
      pwMat[nms[pw[2,i]], nms[pw[1,i]]] <- ovrallpw[i]
    }
  }
  # write results
  if(!is.null(outfile)){
    old_scipen <- options()$scipen
    old_opt <- options(scipen = old_scipen) 
    on.exit(options(old_opt))
    options(scipen=999)
    opf <- paste(getwd(), "/", outfile, "-[chiCalc]/", sep = "")
    dir.create(opf, showWarnings = FALSE)
    glbOut <- c("Loci\tp.value", apply(glb, 1, paste, collapse = "\t"))
    fh <- paste("Fisher's Exact tests for sample independence (global)",
                "Data used: Genotype counts",
                paste("Monte Carlo replications: ", mcRep, sep = ""),
                "",
                "Locus P Values calculated across all samples.",
                "Overall P Value calculated using Fisher's method.",
                "",
                "",
                sep = "\n")
    glbOut <- c(fh, glbOut)
    glbOut <- paste(glbOut, collapse = "\n")
    writeLines(glbOut, paste(opf, "global_chicalc.txt"))
    rm(glbOut)
    if(pairwise){
      pwMatOut <- format(pwMat, nsmall = 4)
      pwMatOut <- rbind(nms, pwMatOut)
      pwMatOut <- cbind(c("pops", nms), pwMatOut)
      pwMatOut <- gsub("\\s+NA", "", pwMatOut)
      diag(pwMatOut)[-1] <- "--"
      pwMatOut <- apply(pwMatOut, 1, paste, collapse = "\t ")
      pwMatOut <- c("\n\nPairwise Matrix of overall P Values:\n", 
                    pwMatOut, "\n", "Locus and overall P Values\n")
      pwMatOut <- paste(pwMatOut, collapse = "\n")
      pwOutFun <- function(pwdat, nms){
        hdr <- paste("\nPairwise comparison: ", nms, "\n", sep = "")
        xOut <- c(hdr, "Loci\tp.value", apply(pwdat, 1, paste, 
                                              collapse = "\t"))
        return(paste(xOut, collapse = "\n"))
      }
      pwOut <- mapply(pwOutFun, pwdat = pwDif, nms = pwnms, SIMPLIFY = FALSE)
      fh <- paste("Fisher's Exact tests for sample independence (Pairwise)",
                  "Data used: Genotype counts",
                  paste("Monte Carlo replications: ", mcRep, sep = ""),
                  "",
                  "Locus P Values calculated across all samples.",
                  "Overall P Value calculated using Fisher's method.",
                  "",
                  "",
                  pwMatOut,
                  sep = "\n")
      pwOut <- paste(c(fh, do.call("c", pwOut)), collapse = "\n")
      writeLines(pwOut, paste(opf, "pairwise_chicalc.txt"))
    }
  }
  
  if(pairwise){
    # generate output
    pwglb <- data.frame(pops = pwnms, p.value = ovrallpw)
    pwDif <- lapply(pwDif, function(x){
      x$p.value[-(which(x$loci == "Overall"))]
    })
    pwDif <- as.data.frame(do.call("cbind", pwDif),
                           row.names = dat$locs, 
                           col.names = pwnms)
    list(overall = glb,
         multilocus_pw = pwglb,
         locus_pw = pwDif)
  } else {
    list(overall = glb)
  }
}