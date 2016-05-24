allele1to2 <-
  function(geno, marker.label = NULL, miss.val = NA){
    geno <- as.matrix(geno)
    geno[!(geno %in% c(0,1,2))] <- NA
    if (!all(is.na(miss.val))) {geno[geno %in% miss.val]  <- NA}
    ff <- function(x) {
      a1 <- ifelse(is.na(x), NA, ifelse(x == 0, 1, ifelse(x == 2, 2, 2)))
      a2 <- ifelse(is.na(x), NA, ifelse(x == 0, 1, ifelse(x == 2, 2, 1)))
      cbind(a1, a2)
    }
    gg <- lapply(as.list(as.data.frame(geno)), ff)
    gg <- do.call("cbind", gg)
    gg <- as.matrix(gg)
    if(is.null(marker.label)) {
      marker.label <-  paste("S", 1:(ncol(gg) / 2), sep = "") 
    }
    colnames (gg) <- paste(rep(marker.label, each = 2),
                           rep(c("a1","a2"), ncol(gg) / 2),
                           sep = ".")
    rownames (gg) <- rownames(geno)
    return(gg)
  }
