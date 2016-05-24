
tab.disjonctif.NA <-
  function (tab,rename.level=FALSE) {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind<-(1:n) + n * (unclass(moda) - 1)
      indNA<-which(is.na(ind))
      
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA,]<-NA 
      if (rename.level==TRUE){
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), 
                                                  sep = "="))
      }    else{
        dimnames(x) <- list(row.names(tab), levels(moda))
      }
      
      
      return(x)
    }
    if (ncol(tab) == 1) 
      res <- modalite.disjonctif(1)
    else {
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
  }

