num.to.lett <- function (xx) 
{
  s1 <- seq(1, 5000, by = 2)
  s2 <- seq(2, 5000, by = 2)
  mark.list <- list(NA)
  names.list <- list(NA)
  for (i in 1:(dim(xx)[2]/2)) {
    v1 <- s1[i]
    v2 <- s2[i]
    names.list[[i]] <- names(xx)[v1]
    xxx <- xx[, v1:v2]
    wer <- unique(which(is.na(xxx), arr.ind=T)[,1])
    if(length(wer) > 0){
      geno <- table(paste(xxx[-wer, 1], xxx[-wer, 2], sep = "-"))/sum(table(paste(xxx[-wer, 
                                                                                      1], xxx[-wer, 2], sep = "-")))
    }else{
      geno <- table(paste(xxx[, 1], xxx[, 2], sep = "-"))/sum(table(paste(xxx[, 
                                                                              1], xxx[, 2], sep = "-")))
    }
    geno <- sort(geno, decreasing = T)
    v <- which(geno < 0.1)
    if (length(v) > 0) {
      geno <- geno[-v]
    }
    alle <- table(c(xxx[, 1], xxx[, 2], sep = "-"))/sum(table(c(xxx[, 
                                                                    1], xxx[, 2], sep = "-")))
    alle <- sort(alle, decreasing = T)
    v <- which(alle < 0.1)
    if (length(v) > 0) {
      alle <- alle[-v]
    }
    n.all <- length(names(alle))
    alls <- as.numeric(names(alle))
    config <- paste(xxx[, 1], xxx[, 2], sep = "-")
    y <- vector(mode = "character", length = length(config))
    if (length(alls) == 2 & length(geno) == 2) {
      y[which(config == paste(alls[1], alls[2], sep = "-") | 
                config == paste(alls[2], alls[1], sep = "-"))] <- "R"
      y[which(config == paste(alls[1], alls[1], sep = "-"))] <- "A"
      y[which(config == paste(alls[2], alls[2], sep = "-"))] <- "A"
    }
    if (length(alls) == 2 & length(geno) == 3) {
      y[which(config == paste(alls[1], alls[1], sep = "-"))] <- "A"
      y[which(config == paste(alls[2], alls[1], sep = "-") | 
                config == paste(alls[1], alls[2], sep = "-"))] <- "R"
      y[which(config == paste(alls[2], alls[2], sep = "-"))] <- "G"
    }
    if (length(alls) == 3 & length(geno) == 4) {
      y[which(config == paste(alls[1], alls[1], sep = "-"))] <- "A"
      y[which(config == paste(alls[2], alls[1], sep = "-") | 
                config == paste(alls[1], alls[2], sep = "-"))] <- "R"
      y[which(config == paste(alls[3], alls[1], sep = "-") | 
                config == paste(alls[1], alls[3], sep = "-"))] <- "W"
      y[which(config == paste(alls[3], alls[2], sep = "-") | 
                config == paste(alls[2], alls[3], sep = "-"))] <- "K"
    }
    if (length(alls) == 4) {
      y[which(config == paste(alls[1], alls[1], sep = "-"))] <- "A"
      y[which(config == paste(alls[2], alls[2], sep = "-"))] <- "T"
      y[which(config == paste(alls[3], alls[3], sep = "-"))] <- "C"
      y[which(config == paste(alls[4], alls[4], sep = "-"))] <- "G"
      y[which(config == paste(alls[2], alls[1], sep = "-") | 
                config == paste(alls[1], alls[2], sep = "-"))] <- "W"
      y[which(config == paste(alls[3], alls[1], sep = "-") | 
                config == paste(alls[1], alls[3], sep = "-"))] <- "M"
      y[which(config == paste(alls[4], alls[1], sep = "-") | 
                config == paste(alls[1], alls[4], sep = "-"))] <- "R"
      y[which(config == paste(alls[2], alls[3], sep = "-") | 
                config == paste(alls[3], alls[2], sep = "-"))] <- "Y"
      y[which(config == paste(alls[2], alls[4], sep = "-") | 
                config == paste(alls[4], alls[2], sep = "-"))] <- "K"
      y[which(config == paste(alls[3], alls[4], sep = "-") | 
                config == paste(alls[3], alls[4], sep = "-"))] <- "S"
    }
    y[which(y == "")] <- NA
    mark.list[[i]] <- y
  }
  xx2 <- data.frame(matrix(unlist(mark.list), nrow = dim(xx)[1]))
  rownames(xx2) <- rownames(xx)
  names(xx2) <- unlist(names.list)
  return(xx2)
}