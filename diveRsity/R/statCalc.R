#' ### StatCalc: Function for the calculation of allele frequencies etc.
#' 
#' __Kevin Keenan__ (2014)
#' 
#' Subfunction for diffCalcRcpp
#' # Calculate pre-diversity/differentiation statistics
#' 
#' This function accepts the output from the `rpg` function and returns
#' some basic statistics which can then be used to calculate Fst etc.
#' 
#' If `statCalc` is passed `idx`, the index rows will be sampled from `plAr`
#' and the statistics describing this re-sample will be returned. Thus, the 
#' function can be used for bootstrapping data.
#' 
# myTab <- function(x){
#   x <- as.character(na.omit(x))
#   mtch <- unique(x)
#   return(sapply(mtch, function(y){
#     return(sum(x == y))
#   }))
# }

# library(Rcpp)
# sourceCpp("src/myTab.cpp")
#' __Kevin Keenan__ (2014)
####-- prestats function --####
statCalc <- function(rsDat, idx = NULL, al, fst, bs = TRUE, 
                     ci_type = "individuals"){
  # generate resamples
  if(bs & ci_type == "individuals"){
    rsFun <- function(x, y){
      return(x[y,,])
    }
    rsDat <- mapply(rsFun, x = rsDat, y = idx, SIMPLIFY = FALSE) 
  } else if(bs & ci_type == "loci"){
    rsFun <- function(x, y){
      return(x[,y,])
    }
    rsDat <- lapply(rsDat, rsFun, y = idx)
    al <- al[idx]
  }
  
  # calculate allele frequecies
  
  alf <- lapply(rsDat, function(x){
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
  nloci <- length(al)
  
  # organise allele frequencies
  alf <- lapply(1:length(al), function(i){
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
  alOut <- mapply(alSort, x = alf, y = al, SIMPLIFY = FALSE)
  # calculate harmonic sizes (only for estimators)
  popSizes <- lapply(rsDat, function(x){
    lgths <- apply(x, 2, function(y){
      nrow(na.omit(y))
    })
    return(lgths)
  })
  ps <- do.call("cbind", popSizes)
  
  # inds per locus
  indtyp <- lapply(1:nloci, function(i){
    vapply(rsDat, function(x){
      op <- length(x[!is.na(x[,i,1]),i,1])
      if(op == 0L){
        return(NA)
      } else {
        return(op) 
      }
    }, FUN.VALUE = numeric(1))
  })
  
  # if wc fst = true calculate
  if(fst){
    hsums <- lapply(rsDat, function(x){
      hts <- lapply(1:dim(x)[2], function(i){
        out <- x[,i,1] == x[,i,2]
        return(out)
      })
      gts <- lapply(1:dim(x)[2], function(i){x[,i,]})
      alls <- lapply(gts, function(y){unique(y[!is.na(y)])})
      #     htcount <- function(gts, hts){
      #      return(table(gts[!hts,]))
      #     }
      htCount <- function(gts, hts, alls){
        if(all(is.na(gts))){
          out <- 0
          names(out) <- "NA"
          return(out)
        } else {
          gts <- gts[!hts, ]
          ht <- sapply(alls, function(al){
            return(sum((gts == al), na.rm = TRUE))
          })
          names(ht) <- alls
          return(ht)
        }
      }
      htcounts <- mapply(htCount, gts = gts, hts = hts, alls = alls, 
                         SIMPLIFY = FALSE)
      return(htcounts)
    })
    # convert hsums to locus focus
    hsums <- lapply(1:nloci, function(i){
      lapply(hsums, "[[", i)
    })
    list(alOut = alOut, ps = ps, hsums = hsums, indtyp = indtyp)
  } else {
    list(alOut = alOut, ps = ps, indtyp = indtyp)
  }
}
##########
# END
##########