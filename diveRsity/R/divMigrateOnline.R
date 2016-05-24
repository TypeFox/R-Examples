# divMigrateOnline function definition for online application
# function definition ----
divMigrateOnline <- function(infile = NULL, nbs = 0, stat = "all", 
                             para = FALSE, true.nm = FALSE, alpha = 0.05){
  # preabmle ----
  #nbs <- 10
  #cat("The method used in this function is still under development. \n")
  # read data ----
  #source("bsFun-fix.R")
  #Rcpp::sourceCpp("pwHt-fix.cpp")
  #bsFun <- diveRsity:::bsFun
  #pwHt <- diveRsity:::pwHt
  #rgp <- diveRsity:::rgp
  #nbs <- 500
  #stat = "all"
  #filter_threshold <- 0
  #plot_network = TRUE
  #plot_col <- "darkblue"
  #para = FALSE
  #infile <- "Inoue_et_al_2013.gen"#"YOSEonly.gen"
  #outfile <- NULL
  #data(Test_data, package = "diveRsity")
  #Test_data[is.na(Test_data)] <- ""
  #Test_data[Test_data == "0"] <- "000000"
  dat <- rgp(infile)
  npops <- length(dat$genos)
  nloci <- length(dat$af)
  # fix allele frequencies
  dat$af <- lapply(dat$af, function(x){
    cs <- colSums(x)
    x[,cs == 0] <- NA
    return(x)
  })
  # generate pw combos ----
  pw <- combn(npops, 2)
  # extract population names (first individual of each sample)
  popnms <- sapply(dat$indnms, "[", 1)
  # calculate ht and hs ----
  #library(Rcpp) # comment out for package
  #sourceCpp("src/pwHt.cpp") # comment out for package
  hths <- lapply(dat$af, pwHt, pw = pw-1)
  # seperate ht and hs matrices
  ht <- lapply(hths, "[[", "ht")
  hs <- lapply(hths, "[[", "hs")
  # replace NaN with NA
  #ht <- lapply(ht, function(x){ x[is.nan(x)] <- NA; return(x)})
  #hs <- lapply(hs, function(x){ x[is.nan(x)] <- NA; return(x)})
  # Calculate D ----
  # function for locus d
  #if(stat == "d" || stat == "Nm"){
    d <- function(ht, hs){
      return(((ht-hs)/(1-hs))*2)
    }
  #}
  # Gst function
  #if(stat == "gst" || stat == "Nm"){
    g <- function(ht, hs){
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
  #}
  # Nm estimator (Alcala et al 2014)
  #if(stat == "Nm"){
    Nm <- function(g, d, n){
      t1 <- (1-g)/g
      t2 <- ((n-1)/n)^2
      t3 <- ((1-d)/(1-((n-1)/n)*d))
      return(0.25*t1*t2*t3)
    }
  #}
  
  # D calculations ----
  #if(stat == "d" || stat == "Nm"){
    dloc <- mapply(`d`, ht = ht, hs = hs, SIMPLIFY = "array")
    dloc[is.nan(dloc)] <- 1
    hrmD <- apply(dloc, c(1,2), function(x){
      mn <- mean(x, na.rm = TRUE)
      vr <- var(x, na.rm = TRUE)
      return(1/((1/mn) + vr * (1/mn)^3))
    })
    dMig <- (1 - hrmD) / hrmD
    # fix infinities
    dMig[is.infinite(dMig)] <- NA
    # calculate relative migration
    dRel <- dMig/max(dMig, na.rm = TRUE)
    dRel[is.nan(dRel)] <- NA
    colnames(dRel) <- popnms
    rownames(dRel) <- popnms
  #}
  
  # Gst calculations ----
  #if(stat == "gst" || stat == "Nm"){
    g <- function(ht, hs){
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
    hsAr <- array(unlist(hs), dim = c(npops, npops, nloci))
    mnHs <- apply(hsAr, c(1,2), mean, na.rm = TRUE)
    htAr <- array(unlist(ht), dim = c(npops, npops, nloci))
    mnHt <- apply(htAr, c(1,2), mean, na.rm = TRUE)
    hrmGst <- g(mnHt, mnHs)
    # calculate migrations from Gst
    gMig <- ((1/hrmGst) - 1)/4
    gMig[is.infinite(gMig)] <- NA
    gRel <- gMig/max(gMig, na.rm = TRUE)
    colnames(gRel) <- popnms
    rownames(gRel) <- popnms
  #}
  
  # Nm calculations ----
  #if(stat == "Nm"){
    nm <- Nm(hrmGst, hrmD, 2)
    diag(nm) <- NA
    nmRel <- nm/max(nm, na.rm = TRUE)
    colnames(nmRel) <- popnms
    rownames(nmRel) <- popnms
    if(true.nm){
      colnames(nm) <- popnms
      rownames(nm) <- popnms
    }
  #}
  
  # Bootstrapping ----
  if(nbs != 0L){
    # generate bootstrap indexes ----
    ps <- sapply(dat$indnms, length)
    idx <- lapply(1:nbs, function(i){
      lapply(ps, function(x){
        return(sample(x, size = x, replace = TRUE))
      })
    })
    
    # calculate bootstrap D ----
    # load bs function
    #source("R/bsFun.R")
    # run bootstrap function
    if(para){
      
      cl <- parallel::makeCluster(detectCores())
      parallel::clusterExport(cl, c("bsFun", "dat", "pw", "stat"),
                              envir = environment())
      bsStat <- parallel::parLapply(cl, idx, function(x){
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, pw = pw,
                     stat = stat))
      })
      parallel::stopCluster(cl)
    } else {
      bsStat <- lapply(idx, function(x){
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, pw = pw,
                     stat = stat))
      })
    }
    
    # convert stats to arrays
    #if(stat == "d"){
      bsD <- sapply(bsStat, "[[", "dRel", simplify = "array")
    #}
    #if(stat == "gst"){
      bsG <- sapply(bsStat, "[[", "gRel", simplify = "array")
    #}
    #if(stat == "Nm"){
      bsNm <- sapply(bsStat, "[[", "nmRel", simplify = "array")
    #}
    
    # function for significant difference determination
    sigDiff <- function(x, y){
      if(x[1] < y[1] && x[2] < y[1]){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    #if(stat == "d"){
      sigMatD <- matrix(NA, nrow = ncol(dRel), ncol(dRel))
      for(i in 1:ncol(pw)){
        p1 <- quantile(bsD[pw[1,i], pw[2,i],], prob = c(alpha/2, 1-(alpha/2)))
        p2 <- quantile(bsD[pw[2,i], pw[1,i],], prob = c(alpha/2, 1-(alpha/2)))
        sigMatD[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
        sigMatD[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
      }
      dRelbs <- dRel
      dRelbs[!sigMatD] <- 0
    #}
    #if(stat == "gst"){
      sigMatG <- matrix(NA, nrow = ncol(gRel), ncol(gRel))
      for(i in 1:ncol(pw)){
        p1 <- quantile(bsG[pw[1,i], pw[2,i],], prob = c(alpha/2, 1-(alpha/2)))
        p2 <- quantile(bsG[pw[2,i], pw[1,i],], prob = c(alpha/2, 1-(alpha/2)))
        sigMatG[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
        sigMatG[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
      }
      gRelbs <- gRel
      gRelbs[!sigMatG] <- 0
    #}
    #if(stat == "Nm"){
      sigMatNm <- matrix(NA, nrow = ncol(nmRel), ncol(nmRel))
      for(i in 1:ncol(pw)){
        p1 <- quantile(bsNm[pw[1,i], pw[2,i],], prob = c(alpha/2, 1-(alpha/2)))
        p2 <- quantile(bsNm[pw[2,i], pw[1,i],], prob = c(alpha/2, 1-(alpha/2)))
        sigMatNm[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
        sigMatNm[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
      }
      nmRelbs <- nmRel
      nmRelbs[!sigMatNm] <- 0
    #}
  }
    if(nbs != 0L && !true.nm){
      list(D = round(dRel, 3),
           D_sig = sigMatD,
           Gst = round(gRel, 3),
           G_sig = sigMatG,
           Nm = round(nmRel, 3),
           Nm_sig = sigMatNm,
           nbs = nbs)
    } else if(nbs != 0L && true.nm){
      list(D = round(dRel, 3),
           D_sig = sigMatD,
           Gst = round(gRel, 3),
           G_sig = sigMatG,
           Nm = round(nmRel, 3),
           Nm_sig = sigMatNm,
           true_nm = nm,
           nbs = nbs)
    } else if(nbs == 0L && !true.nm){
      list(D = round(dRel, 3),
           Gst = round(gRel, 3),
           Nm = round(nmRel, 3),
           nbs = 0L)
    } else if(nbs == 0L && true.nm){
      list(D = round(dRel, 3),
           Gst = round(gRel, 3),
           Nm = round(nmRel, 3),
           true_nm = nm,
           nbs = 0L)
    }
}