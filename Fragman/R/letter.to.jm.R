letter.to.jm <- function (x) 
{
  allelesf <- na.omit(unique(x))
  all1 <- table(na.omit((x)))/length(x)
  v <- which(all1 > 0.1)
  allelesff <- names(all1)[v]
  v2 <- which(as.vector(unlist(allelesf)) %in% allelesff)
  alleles <- as.vector(unlist(allelesf))[v2]
  if (length(alleles) == 1) {
    alleles1 <- sort(alleles)
    y <- rep(NA, length(x))
    y[which(x == alleles1[1])] <- "-"
    y[which(is.na(x))] <- "--"
  }
  if (length(alleles) == 2) {
    alleles1 <- alleles
    if ((alleles1[1] == "A") | (alleles1[1] == "C") | (alleles1[1] == 
                                                       "T") | (alleles1[1] == "G")) {
      y <- rep(NA, length(x))
      y[which(x == alleles1[1])] <- "nn"
      y[which(x == alleles1[2])] <- "np"
      y[which(is.na(x))] <- "--"
    }
    if ((alleles1[1] == "R") | (alleles1[1] == "Y") | (alleles1[1] == 
                                                       "S") | (alleles1[1] == "W") | (alleles1[1] == "K") | 
        (alleles1[1] == "M")) {
      y <- rep(NA, length(x))
      y[which(x == alleles1[1])] <- "lm"
      y[which(x == alleles1[2])] <- "ll"
      y[which(is.na(x))] <- "--"
    }
  }
  if (length(alleles) == 3) {
    yuyu <- c(length(which(x == alleles[1])), length(which(x == 
                                                             alleles[2])), length(which(x == alleles[3])))
    names(yuyu) <- alleles
    alleles1 <- names(sort(yuyu, decreasing = T))
    y <- rep(NA, length(x))
    y[which(x == alleles1[1])] <- "hk"
    y[which(x == alleles1[2])] <- "hh"
    y[which(x == alleles1[3])] <- "kk"
    y[which(is.na(x))] <- "--"
  }
  if (length(alleles) == 4) {
    alleles1 <- alleles
    # R = A:G
    # K = G:T
    # A = A:A
    # W = A:T
    #alleles2 <- sort(alleles)
    mat <- alleles1[which(alleles1 %in% x[1])]
    pat <- alleles1[which(alleles1 %in% x[2])]
    
    #oc <- c(alleles1[1], alleles1[2], alleles2[1])
    #v <- which(alleles1 == oc)
    #v2 <- 1:4
    #v3 <- which(v2 != v)
    y <- rep(NA, length(x))
    y[which(x == mat)] <- "ef"
    y[which(x == pat)] <- "eg"
    y[which(x == "A")] <- "ee"
    y[which(x == "K")] <- "fg"
    y[which(is.na(x))] <- "--"
  }
  if (length(alleles) == 6) {
    alleles1 <- alleles
    if (alleles1[1] == "R") {
      y <- rep(NA, length(x))
      y[which(x == "M")] <- "ac"
      y[which(x == "W")] <- "ad"
      y[which(x == "S")] <- "bc"
      y[which(x == "K")] <- "bd"
      y[which(x == "R")] <- "ab"
      y[which(x == "Y")] <- "cd"
      y[which(is.na(x))] <- "--"
    }
    if (alleles1[1] == "Y") {
      y <- rep(NA, length(x))
      y[which(x == "M")] <- "ac"
      y[which(x == "W")] <- "bc"
      y[which(x == "S")] <- "ad"
      y[which(x == "K")] <- "bd"
      y[which(x == "Y")] <- "ab"
      y[which(x == "R")] <- "cd"
      y[which(is.na(x))] <- "--"
    }
    if (alleles1[1] == "S") {
      y <- rep(NA, length(x))
      y[which(x == "M")] <- "ac"
      y[which(x == "Y")] <- "ad"
      y[which(x == "R")] <- "bc"
      y[which(x == "K")] <- "bd"
      y[which(x == "S")] <- "ab"
      y[which(x == "W")] <- "cd"
      y[which(is.na(x))] <- "--"
    }
    if (alleles1[1] == "W") {
      y <- rep(NA, length(x))
      y[which(x == "M")] <- "ac"
      y[which(x == "Y")] <- "bc"
      y[which(x == "R")] <- "ad"
      y[which(x == "K")] <- "bd"
      y[which(x == "W")] <- "ab"
      y[which(x == "S")] <- "cd"
      y[which(is.na(x))] <- "--"
    }
    if (alleles1[1] == "K") {
      y <- rep(NA, length(x))
      y[which(x == "R")] <- "ac"
      y[which(x == "S")] <- "ad"
      y[which(x == "W")] <- "bc"
      y[which(x == "Y")] <- "bd"
      y[which(x == "K")] <- "ab"
      y[which(x == "M")] <- "cd"
      y[which(is.na(x))] <- "--"
    }
    if (alleles1[1] == "M") {
      y <- rep(NA, length(x))
      y[which(x == "R")] <- "ac"
      y[which(x == "S")] <- "bc"
      y[which(x == "W")] <- "ad"
      y[which(x == "Y")] <- "bd"
      y[which(x == "M")] <- "ab"
      y[which(x == "K")] <- "cd"
      y[which(is.na(x))] <- "--"
    }
  }
  y[which(is.na(y) == TRUE)] <- "--"
  return(y)
}