#' divInd function for calculating Kosman & Leonard (2005) similarity 
#' coefficients in diploid individuals.
#' 
#' Kevin Keenan, 2014

#' divInd function for calculating Kosman & Leonard (2005) similarity 
#' coefficients in diploid individuals.
#' 
#' Kevin Keenan, 2014

divSimCo <- function(infile = NULL, loci = FALSE, boots = 0){
  #infile <- "monomorphism_main.gen"
  #boots = 10
  #loci = FALSE
  #genos2mat <- diveRsity::genos2mat
  dat <- rgp(infile)
  nms <- gsub(",", "", unlist(dat$indnms))
  stats <- sapply(dat$genos, function(x) dim(x)[1:2])
  
  dims <- c(rowSums(stats)[1], stats[2,1])
  
  genos <- array(NA, dim = c(dims, 2))
  for(i in 1:ncol(stats)){
    if(i == 1L){
      j = 1
    }
    idx <- sum(stats[1,(1:i)])
    genos[j:idx,,] <- dat$genos[[i]]
    j = idx+1
  }
  
  # get all alleles per locus
  test <- apply(genos, 2, function(x){
    nms <- as.character(na.omit(unique(unlist(apply(x, 2, unique)))))
    inmat <- array(0, dim = c(nrow(x), length(nms)))
    na <- rep(NA, length(nms))
    Reduce(`+`, lapply(1:2, function(i){
      ip <- match(x[,i], nms)-1
      return(genos2mat(mat = inmat, ip, na))
    }))
  })
  
  # calculate the distance for each individual
  if(loci || !boots == 0L){
    
    dists <- lapply(test, function(x){
      as.matrix(dist(x, method = "manhattan"))/2
    })
    
    dists2 <- lapply(dists, function(x){
      x[is.na(x)] <- 0
      return(x)
    })
    
    rm(test)
    
    indLoc <- colSums(apply(genos[,,1], 1, function(x) !is.na(x)))
    
    
    glb <- Reduce(`+`, dists2) %*% diag(1/indLoc)
    
    rm(dists2)
  } else {
    dists <- lapply(test, function(x){
      out <- as.matrix(dist(x, method = "manhattan"))/2
      out[is.na(out)] <- 0
      return(out)
    })
    
    indLoc <- colSums(apply(genos[,,1], 1, function(x) !is.na(x)))
    
    glb <- Reduce(`+`, dists) %*% diag(1/indLoc)
    
    rm(dists)
  }
  
  # bootstrapping code. Not finalised!
  if(boots > 0L){
    
    idx <- lapply(1:boots, function(i){
      sample(dim(genos)[2], dim(genos)[2], replace = TRUE)
    })
    
    boot_dist <- lapply(idx, function(x){
      loc <- dists[x]
      
      # count the loci scored
      indtyp <- indLoc <- colSums(apply(genos[,x,1], 1, function(x) !is.na(x)))
      
      # replace NAs with 0
      loc <- lapply(loc, function(y){
        y[is.na(y)] <- 0
        return(y)
      })
      
      # calculate multilocus par
      out <- Reduce(`+`, loc) %*% diag(1/indtyp)
      dimnames(out) <- list(nms, nms)
      return(out)
    })
  }
  
  
  dimnames(glb) <- list(nms, nms)
  if(loci){
    dists <- lapply(dists, function(x){
      dimnames(x) <- list(nms, nms)
      return(x)
    })
    names(dists) <- dat$locs
  }
  if(loci & boots == 0L){
    list(glb = glb,
         loci = dists)
  } else if(!loci & boots == 0L){
    list(glb = glb)
  } else if(!loci & !boots == 0L){
    list(glb = glb,
         bootstrapped = boot_dist)
  } else if(loci & !boots == 0L){
    list(glb = glb,
         loci = dists,
         bootstrapped = boot_dist)
  } else {
    list(glb = glb)
  }
}