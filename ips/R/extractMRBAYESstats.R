extractMRBAYESstats <- function(file, nn) {
  
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  
    # isolate NEWICK string
    # ------------------------------
    patt <- "tree.*\\[&U\\]"
    X <- X[grep(patt, X)]
    X <- gsub(patt, "", X)
    X <- gsub(";$", ");", X)
    
    # store stats per node in a list
    # ------------------------------
    tab <- unlist(strsplit(X, "\\["))[-1]
    tab <- gsub("&|;|\\]", "", tab)
    tab <- gsub(":.+$", "", tab)
    tab <- gsub("(,)([[:alpha:]])", "xxx\\2", tab)
    foo <- function(x){x <- unlist(strsplit(x, "xxx")); x}
    tab <- lapply(tab, foo)
    foo2 <- function(x) {
      nms <- gsub("(^)(.*)(=)(.*$)", "\\2", x)
      val <- strsplit(gsub("^.*=", "", x), ",")
      names(val) <- nms
      val <- val[!names(val) %in% c("prob(percent)", "prob+-sd")]
      val <- lapply(val, gsub, pattern = "^[[:punct:]]*|[[:punct:]]*$", 
                    replacement = "")
      lapply(val, as.numeric)
    }
    tab <- lapply(tab, foo2)
    
    tab <- lapply(tab, unlist)
    
    
    id <- seq(1, nn * 2 - 1, by = 2)
    nodestat <- lapply(tab[id], unlist)
    edgestat <- lapply(tab[id + 1], unlist)
    edgestat <- lapply(edgestat, function(x) x[1:4])
    tab <- cbind(do.call(rbind, edgestat), do.call(rbind, nodestat))
    
    return(tab)
}