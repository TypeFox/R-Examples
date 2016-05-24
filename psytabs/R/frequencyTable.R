frequencyTable <-
  function (factors) {
    tab.list <- lapply(factors, freq)
    tab <- data.frame(do.call(rbind, tab.list, 2))
    
    pattern <- unlist(strsplit(rownames(tab)[grepl("X\\[\\[", rownames(tab))], "\\."))
    pattern <- gsub(paste(pattern, collapse=".|"), "", rownames(tab))
    pattern[grepl("X\\[\\[", rownames(tab))] <- names(tab.list)
    rownames(tab) <- pattern
    tab
  }

