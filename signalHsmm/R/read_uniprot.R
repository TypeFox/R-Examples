#' Read data from UniProt database
#'
#' Read data saved in UniProt original flat text format.
#'
#' @param connection a \code{\link{connection}} to UniProt data in text format.
#' @param ft_names a character vector of UuniProt features to be extracted, for example 
#' \code{"signal"}, \code{"transit"}, \code{"propep"}. The case is not matched.
#' @param kwds a \code{NULL} or character vector of keywords (not UniProt keywords, 
#' but words of interest, that may occur in the protein description).
#' @keywords manip
#' @return a list of sequences. Each element has a class \code{\link[seqinr]{SeqFastaAA}}.
#' Attributes \code{OS} and \code{OC} represents respectively OS and OC fields in the 
#' protein description. A value of each feature is preserved as an attribute named 
#' after the feature.
#' @export
#' @keywords manip

read_uniprot <- function(connection, ft_names, kwds = NULL) {
  
  all_lines <- readLines(connection)
  
  prot_ids <- grep("\\<ID   ", all_lines)
  
  # get features
  fts <- lapply(ft_names, function(single_ft) {
    get_ft(all_lines, prot_ids, single_ft)
  })
  
  # validate features
  validated_prots <- apply(sapply(fts, validate_ft, donts_symbols = c(">", "<1", "?", "Or ", "Not cleaved")), 
                           1, all)
  #validation results
  valres <- matrix(c(prot_ids, validated_prots), ncol = 2, dimnames = list(NULL, c("id", "selected")))
  
  if(!is.null(kwds)) {
    kwd_list <- do.call(cbind, lapply(kwds, function(single_ft) {
      get_keyword(all_lines, prot_ids, single_ft)
    }))
    valres[, "selected"] <- valres[, "selected"] & apply(kwd_list, 1, any)
  }
  
  all_seqs <- cbind(valres, get_add_id(all_lines))
  
  sure_seqs <- all_seqs[all_seqs[, "selected"] == 1, ]
  
  list_prots <- lapply(1L:nrow(sure_seqs), function(i) {
    start_seq <- sure_seqs[i, "seq_start"]
    end_seq <- sure_seqs[i, "seq_end"]
    
    ith_seq <- strsplit(gsub(" ", "", paste0(all_lines[start_seq:end_seq], collapse = "")), "")[[1]]
    
    class(ith_seq) <- "SeqFastaAA"
    aa_name <- strsplit(all_lines[sure_seqs[i, "id"]], "   ")[[1]][2]
    attr(ith_seq, "name") <- aa_name
    attr(ith_seq, "Annot") <- paste0(">", aa_name)
    attr(ith_seq, "class") <- "SeqFastaAA"
    
    os <- all_lines[sure_seqs[i, "os_start"]:sure_seqs[i, "os_end"]]
    attr(ith_seq, "OS") <- paste0(vapply(strsplit(os, "OS   "), function(single_os) 
      single_os[2], "a"), collapse = " ")
    
    oc <- all_lines[sure_seqs[i, "oc_start"]:sure_seqs[i, "oc_end"]]
    attr(ith_seq, "OC") <- paste0(vapply(strsplit(oc, "OC   "), function(single_oc) 
      single_oc[2], "a"), collapse = " ")
    
    ith_seq
  })
  
  for(feature_id in 1L:length(ft_names)) {
    #use only validated features
    sure_fts <- fts[[feature_id]][all_seqs[, "selected"] == 1]
    for(i in 1L:length(list_prots)) {
      attr(list_prots[[i]], ft_names[feature_id]) <- as.numeric(strsplit(sure_fts[[i]], "[ ]+")[[1]][3L:4])
    }
  }
  
  names(list_prots) <- sapply(list_prots, function(i) attr(i, "name"))
  atypical <- get_atyp(list_prots)
  if (length(atypical) > 0){
    list_prots[-atypical]
  } else {
    list_prots
  }
}

#READ UNIPROT DATA -----------------------------
#helper function to get seqs from  .txt files
#get ids additional data: sequence of amino acids, os, oc field
get_add_id <- function(all_lines) {
  
  seqs_start <- grep("\\<SQ   ", all_lines) + 1
  seqs_end <-  grep("^//", all_lines) - 1
  
  os_list <- get_block(grep("\\<OS   ", all_lines))
  os_start <- vapply(os_list, function(i) i[1], 0)
  os_end <- vapply(os_list, function(i) i[length(i)], 0)
  
  
  oc_list <- get_block(grep("\\<OC   ", all_lines))
  oc_start <- vapply(oc_list, function(i) i[1], 0)
  oc_end <- vapply(oc_list, function(i) i[length(i)], 0)
  
  matrix(c(seqs_start, seqs_end, os_start, os_end, oc_start, oc_end),
         ncol = 6, dimnames = list(NULL, c("seq_start", "seq_end", 
                                           "os_start", "os_end",
                                           "oc_start", "oc_end")))
  
}

#get lines with feature
get_ft <- function(all_lines, prot_ids, ft) {
  #not smart, duplicating long prot_ids vector
  prot_ids <- c(prot_ids, length(all_lines) + 1)
  lapply(2L:length(prot_ids), function(id) {
    sublines <- all_lines[prot_ids[id - 1]:(prot_ids[id] - 1)]
    sublines[grep(paste0("\\<FT   ", toupper(ft)), sublines)]
  })
}

#validate feature
#returns TRUE if feature is positively validated
validate_ft <- function(ft_list, donts_symbols, ft_length = 2) {
  #single feature line
  #donts_symbols symbols which should be not used
  #ft_length feature must be shorter than
  sapply(ft_list, function(single_ftline) {
    (length(single_ftline) < ft_length && !any(vapply(donts_symbols, function(single_symbol)
      grepl(single_symbol, single_ftline, fixed = TRUE), TRUE)))
  })
}



#remove seqs with atypical and/or not identified aas
get_atyp <- function(list_prots) {
  which(colSums(sapply(list_prots, function(prot) 
    sapply(c("B", "U", "X", "Z"), function(atyp_aa)
      atyp_aa %in% prot))) > 0)
}

#gets list of starts and ends of specific blocks in UniProt data
get_block <- function(block_indices) {
  block_list <- list()
  tmp_block <- c()
  
  for(i in block_indices) {
    if(is.null(tmp_block[length(tmp_block)])) {
      tmp_block <- i
    } else {
      if(tmp_block[length(tmp_block)] + 1 == i) {
        tmp_block <- c(tmp_block, i)
      } else {
        block_list <- c(block_list, list(tmp_block))
        tmp_block <- i
      } 
    }
  }
  
  c(block_list, list(tmp_block))
}

#get keyword
get_keyword <- function(all_lines, prot_ids, keyword) {
  prot_ids <- c(prot_ids, length(all_lines) + 1)
  sapply(2L:length(prot_ids), function(id) {
    sublines <- all_lines[prot_ids[id - 1]:(prot_ids[id] - 1)]
    sum(grepl(keyword, sublines)) > 0
  })
}