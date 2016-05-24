#' divMigrate-Fix
#' 
#' Kevin Keenan
#' 
# preamble ----

#' divMigrate: an experimental function for detecting directional differentiation 
#' A function to calculate pairwise directional differentiation
#' a presented in the paper 'Directional genetic differentiation and
#' asymmetric migration Lisa Sundqvist, Martin Zackrisson & David Kleinhans,
#' 2013, arXiv pre-print (http://arxiv.org/abs/1304.0118)'
#' @export

# function definition ----
divMigrate <- function(infile = NULL, outfile = NULL, boots = 0, stat = "all",
                       filter_threshold = 0, plot_network = FALSE, 
                       plot_col = "darkblue", para = FALSE){
  # preabmle ----
  #boots <- 10
  #cat("The method used in this function is still under development. \n")
  # read data ----
  #source("bsFun-fix.R")
  #Rcpp::sourceCpp("pwHt-fix.cpp")
  #pwHt <- diveRsity:::pwHt
  #rgp <- diveRsity:::rgp
  #boots <- 0
  #stat = "all"
  #filter_threshold <- 0
  #plot_network = FALSE
  #plot_col <- "darkblue"
  #para = FALSE
  #data(Test_data, package = "diveRsity")
  #infile <- Test_data
  #outfile = NULL
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
  
  if(!is.null(outfile)){
    dir.create(path = paste(getwd(), "/", outfile, "-[divMigrate]", "/", 
                            sep = ""))
    of <- paste(getwd(), "/", outfile, "-[divMigrate]", "/", sep = "")
  }
  
  # generate pw combos ----
  pw <- combn(npops, 2)
  
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
  if(stat == "d" || stat == "all" || stat == "Nm"){
    d <- function(ht, hs){
      return(((ht-hs)/(1-hs))*2)
    } 
  }
  # Gst function
  if(stat == "gst" || stat == "all" || stat == "Nm"){
    g <- function(ht, hs){
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
  }
  # Nm estimator (Alcala et al 2014)
  if(stat == "Nm" || stat == "all"){
    Nm <- function(g, d, n){
      t1 <- (1-g)/g
      t2 <- ((n-1)/n)^2
      t3 <- ((1-d)/(1-((n-1)/n)*d))
      return(0.25*t1*t2*t3)
    }
  }
  
  # D calculations ----
  if(stat == "d" || stat == "all" || stat == "Nm"){
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
  }
  
  # Gst calculations ----
  if(stat == "gst" || stat == "all" || stat == "Nm"){
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
  }
  
  # Nm calculations ----
  if(stat == "all" || stat == "Nm"){
    nm <- Nm(hrmGst, hrmD, 2)
    diag(nm) <- NA
    nmRel <- nm/max(nm, na.rm = TRUE)
  }
  
  if(plot_network || boots != 0L){
    # plotting network
    if(stat == "d" || stat == "all" || stat == "Nm"){
      dRelPlt <- dRel
      dRelPlt[dRelPlt < filter_threshold] <- 0
    }
    if(stat == "gst" || stat == "all" || stat == "Nm"){
      gRelPlt <- gRel
      gRelPlt[gRelPlt < filter_threshold] <- 0
    }
    if(stat == "all" || stat == "Nm"){
      nmRelPlt <- nmRel
      nmRelPlt[nmRelPlt < filter_threshold] <- 0
    }
  }
  if(plot_network){
    if(stat == "d"){
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, 
                     edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; D)", sep = ""))
      if(!is.null(outfile)){
        pdf(paste(of, "Relative_migration.pdf", sep = ""), paper = "a4r")
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; D)", sep = ""))
      }
    }
    if(stat == "gst"){
      qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, 
                     edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Gst)", sep = ""))
      if(!is.null(outfile)){
        pdf(paste(of, "Relative_migration.pdf", sep = ""), paper = "a4r")
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Gst)", sep = ""))
      }
    }
    if(stat == "Nm"){
      qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, 
                     edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Nm)", sep = ""))
      if(!is.null(outfile)){
        pdf(paste(of, "Relative_migration.pdf", sep = ""), paper = "a4r")
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Nm)", sep = ""))
      }
    }
    if(stat == "all"){
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, 
                     edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; D)", sep = ""))
      qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, 
                     edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Gst)", sep = ""))
      qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, 
                     edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Nm)", sep = ""))
      if(!is.null(outfile)){
        pdf(paste(of, "Relative_migration.pdf", sep = ""), paper = "a4r")
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; D)", sep = ""))
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Gst)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2,2,5,5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Nm)", sep = ""))
      }
    }
    if(boots == 0L && !is.null(outfile)){
      dev.off() 
    } 
  }
  # Bootstrapping ----
  
  if(boots != 0L){
    # generate bootstrap indexes ----
    ps <- sapply(dat$indnms, length)
    idx <- lapply(1:boots, function(i){
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
    if(stat == "d" || stat == "all"){
      bsD <- sapply(bsStat, "[[", "dRel", simplify = "array")
    }
    if(stat == "gst" || stat == "all"){
      bsG <- sapply(bsStat, "[[", "gRel", simplify = "array")
    }
    if(stat == "Nm" || stat == "all"){
      bsNm <- sapply(bsStat, "[[", "nmRel", simplify = "array")
    }
    
    # function for significant difference determination
    sigDiff <- function(x, y){
      if(x[1] < y[1] && x[2] < y[1]){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    if(stat == "d" || stat == "all"){
      sigMatD <- matrix(NA, nrow = ncol(dRel), ncol(dRel))
      for(i in 1:ncol(pw)){
        p1 <- quantile(bsD[pw[1,i], pw[2,i],], prob = c(0.025, 0.975))
        p2 <- quantile(bsD[pw[2,i], pw[1,i],], prob = c(0.025, 0.975))
        sigMatD[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
        sigMatD[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
      }
      dRelPlt[!sigMatD] <- 0
    }
    if(stat == "gst" || stat == "all"){
      sigMatG <- matrix(NA, nrow = ncol(gRel), ncol(gRel))
      for(i in 1:ncol(pw)){
        p1 <- quantile(bsG[pw[1,i], pw[2,i],], prob = c(0.025, 0.975))
        p2 <- quantile(bsG[pw[2,i], pw[1,i],], prob = c(0.025, 0.975))
        sigMatG[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
        sigMatG[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
      }
      gRelPlt[!sigMatG] <- 0
    }
    if(stat == "Nm" || stat == "all"){
      sigMatNm <- matrix(NA, nrow = ncol(nmRel), ncol(nmRel))
      for(i in 1:ncol(pw)){
        p1 <- quantile(bsNm[pw[1,i], pw[2,i],], prob = c(0.025, 0.975))
        p2 <- quantile(bsNm[pw[2,i], pw[1,i],], prob = c(0.025, 0.975))
        sigMatNm[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
        sigMatNm[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
      }
      nmRelPlt[!sigMatNm] <- 0
    }
    # plotting
    if(plot_network){
      # plots to file
      if(stat == "d" && !is.null(outfile)){
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; D method)", sep = ""))
        dev.off()
      }
      if(stat == "gst" && !is.null(outfile)){
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Gst method)", sep = ""))
        dev.off()
      }
      if(stat == "Nm" && !is.null(outfile)){
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Nm method)", sep = ""))
        dev.off()
      }
      if(stat == "all" && !is.null(outfile)){
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; D method)", sep = ""))
        
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Nm method)", sep = ""))
        dev.off()
      }
      # standard plots
      if(stat == "d"){
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; D method)", sep = ""))
      }
      if(stat == "gst"){
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Gst method)", sep = ""))
      }
      if(stat == "Nm"){
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Nm method)", sep = ""))
      }
      if(stat == "all"){
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; D method)", sep = ""))
        
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Gst Method)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                       legend = TRUE, posCol = plot_col, 
                       label.color = plot_col,
                       edge.labels = TRUE, curve = 2.5, mar = c(2,2,5,5))
        title(paste("Significant relative migration network \n (", boots, 
                    " bootstraps; Nm method)", sep = ""))
      }
    }
  }
  
  if(boots != 0L){
    if(stat == "d"){
      list(dRelMig = dRel,
           dRelMigSig = dRelPlt)
    } else if(stat == "gst"){
      list(gRelMig = gRel,
           gRelMigSig = gRelPlt)
    } else if(stat == "Nm"){
      list(nmRelMig = nmRel,
           nmRelMigSig = nmRelPlt)
    } else if(stat == "all"){
      list(dRelMig = dRel,
           dRelMigSig = dRelPlt,
           gRelMig = gRel,
           gRelMigSig = gRelPlt,
           nmRelMig = nmRel,
           nmRelMigSig = nmRelPlt)
    }
  } else {
    if(stat == "d"){
      list(dRelMig = dRel)
    } else if(stat == "gst"){
      list(gRelMig = gRel)
    } else if(stat == "Nm"){
      list(nmRelMig = nmRel)
    } else if(stat == "all"){
      list(dRelMig = dRel,
           gRelMig = gRel,
           nmRelMig = nmRel)
    }
  }
}

# code end ----
#' Bootstrapping function for use with divMigrate
#' 
#' Kevin Keenan (2014)

# Bootstrapping function definition ----
bsFun <- function(genos, idx, af, pw, stat){
  nl <- length(af)
  #myTab <- diveRsity:::myTab
  #pwHt <- diveRsity:::pwHt
  # sub-sample genos ----
  sampleFun <- function(input, idx){
    return(input[idx,,])
  }
  # sub-sample genos
  genos <- mapply(sampleFun, input = genos, idx = idx,
                  SIMPLIFY = FALSE)
  # calculate allele frequencies ----
  #sourceCpp("src/myTab.cpp")
  alf <- lapply(genos, function(x){
    apply(x, 2, function(y){
      if(all(is.na(y))){
        return(NA)
      } else {
        y <- as.vector(na.omit(y))
        nms <- unique(y)[order(unique(y))]
        ot <- myTab(y)
        names(ot) <- nms
        return(ot)
      }
    })
  })
  # organise allele frequencies
  alf <- lapply(1:nl, function(i){
    lapply(alf, "[[", i)
  })
  alSort <- function(x, y){
    idx <- lapply(x, function(z){
      match(names(z), rownames(y))
    })
    for(i in 1:length(idx)){
      y[idx[[i]], i] <- x[[i]]
    }
    return(y)
  }
  # generate allele frequency output
  af <- mapply(alSort, x = alf, y = af, SIMPLIFY = FALSE)
  # calculate hths from boostrapped allele frequencies ----
  hths <- lapply(af, pwHt, pw = pw-1)
  # seperate ht and hs matrices
  ht <- lapply(hths, "[[", "ht")
  hs <- lapply(hths, "[[", "hs")
  
  # D calculations ----
  if(stat == "d" || stat == "all" || stat == "Nm"){
    d <- function(ht, hs){
      return(((ht-hs)/(1-hs))*2)
    }
    # locus d
    dloc <- mapply(`d`, ht = ht, hs = hs, SIMPLIFY = "array")
    # calculate the harmonic mean of Locus D ----
    # set any nan values to 1
    dloc[is.nan(dloc)] <- 1
    # calculate multilocus d
    hrmD <- apply(dloc, c(1,2), function(x){
      mn <- mean(x, na.rm = TRUE)
      vr <- var(x, na.rm = TRUE)
      return(1/((1/mn) + vr * (1/mn)^3))
    })
    dMig <- (1 - hrmD) / hrmD
    dMig[is.infinite(dMig)] <- NA
    dRel <- dMig/max(dMig, na.rm = TRUE)
    diag(dRel) <- NA
  }
  npops <- ncol(hs[[1]])
  # Gst calculations ----
  if(stat == "gst" || stat == "all" || stat == "Nm"){
    g <- function(ht, hs){
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
    hsAr <- array(unlist(hs), dim = c(npops, npops, nl))
    mnHs <- apply(hsAr, c(1,2), mean, na.rm = TRUE)
    htAr <- array(unlist(ht), dim = c(npops, npops, nl))
    mnHt <- apply(htAr, c(1,2), mean, na.rm = TRUE)
    hrmGst <- g(mnHt, mnHs)
    # calculate migrations from Gst
    gMig <- ((1/hrmGst) - 1)/4
    gMig[is.infinite(gMig)] <- NA
    gRel <- gMig/max(gMig, na.rm = TRUE)
  }
  
  # Nm calculations ----
  if(stat == "Nm" || stat == "all"){
    Nm <- function(g, d, n){
      t1 <- (1-g)/g
      t2 <- ((n-1)/n)^2
      t3 <- ((1-d)/(1-((n-1)/n)*d))
      return(0.25*t1*t2*t3)
    }
    nm <- Nm(hrmGst, hrmD, 2)
    diag(nm) <- NA
    nmRel <- nm/max(nm, na.rm = TRUE)
  }
  if(stat == "d"){
    list(dRel = dMig/max(dMig, na.rm = T))
  } else if(stat == "gst"){
    list(gRel = gMig/max(gMig, na.rm = T))
  } else if(stat == "Nm"){
    list(nmRel = nm/max(nm, na.rm = T))
  } else if(stat == "all"){
    list(dRel = dMig/max(dMig, na.rm = T), 
         gRel = gMig/max(gMig, na.rm = T), 
         nmRel = nm/max(nm, na.rm = T))
  }
}