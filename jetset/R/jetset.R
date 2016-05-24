#
# jetset.R
#
# Aron C. Eklund 
# 2011-08-18
#

.checkLookup <- function(query, result) {
    noResult <- query[which(is.na(result))]
    if(length(noResult) > 0) {
      all.noResult <- paste(noResult, collapse = ', ')
      warning("these items are not recognized: ", all.noResult, call. = FALSE)
    }
    wh.ambig <- which(sapply(result, length) > 1)
    if(length(wh.ambig) > 0) {
      merged <- sapply(result[wh.ambig], paste, collapse = ',')
      all.ambig <- paste(query[wh.ambig], '={', merged, '}', sep = '', collapse = '; ')
      warning("these items map to multiple Entrez IDs: ", all.ambig, call. = FALSE)
    }
}

jscores <- function(chip, probeset, eg, symbol, alias, ensembl) {
  stopifnot(length(chip) == 1)
  scores <- get(paste('scores.', chip, sep = ''))
  if(missing(probeset)) probeset <- character(0)
  if(!missing(eg)) {  # eg was supplied
    eg <- as.character(eg)  # just in case
  } else {
    eg <- character(0)
  }
  if(!missing(symbol)) {  # symbol was supplied
    symbol.eg.list <- mget(symbol, org.Hs.egSYMBOL2EG, ifnotfound = list(NA))
    .checkLookup(symbol, symbol.eg.list)
    eg <- c(eg, unlist(symbol.eg.list))
  }
  if(!missing(alias)) {  # alias was supplied
    alias.eg.list <- mget(alias, org.Hs.egALIAS2EG, ifnotfound = list(NA))
    .checkLookup(alias, alias.eg.list)
    eg <- c(eg, unlist(alias.eg.list))
  }
  if(!missing(ensembl)) {  # ensembl was supplied
    ensembl.eg.list <- mget(ensembl, org.Hs.egENSEMBL2EG, ifnotfound = list(NA))
    .checkLookup(ensembl, ensembl.eg.list)
    eg <- c(eg, unlist(ensembl.eg.list))
  }
  if(length(eg) + length(probeset) > 0) {
    eg <- na.omit(eg)
    keep <- (scores$EntrezID %in% eg) | (rownames(scores) %in% probeset)
    scores <- scores[keep, ]
  }
  ### calculate robust and overall scores
  p <- ifelse(chip == 'u133x3p', 300, 600)
  scores$robust <- (1 - (1/p)) ^ scores$process
  scores$overall <- with(scores, specificity * coverage * robust)
  scores$symbol <- as.character(rep(NA, nrow(scores)))
  ### look up symbol
  hasEntrezID <- !is.na(scores$EntrezID)
  scores$symbol[hasEntrezID] <- sapply(mget(scores$EntrezID[hasEntrezID], org.Hs.egSYMBOL, ifnotfound = list(NA)), function(x) x[1])
  scores
}


jmap <- function(chip, eg, symbol, alias, ensembl) {
  stopifnot(length(chip) == 1)
#  scores <- get(paste('scores.', chip, sep = ''))
  nArgs <- sum(!missing(eg), !missing(symbol), !missing(alias), !missing(ensembl))
  if(nArgs == 0) stop("either {'eg', 'symbol', 'alias', or 'ensembl'} must be specified")
  if(nArgs > 1) stop("only one of {'eg', 'symbol', 'alias', or 'ensembl'} can be specified")
  if(!missing(eg)) {  # eg was supplied
    nms <- eg
    eg <- as.character(eg)  # just in case
  } else {
    if(!missing(symbol)) {  # symbol was supplied
      nms <- symbol
      eg.list <- mget(symbol, org.Hs.egSYMBOL2EG, ifnotfound = list(NA))
    }
    if(!missing(alias)) {  # alias was supplied
      nms <- alias
      eg.list <- mget(alias, org.Hs.egALIAS2EG, ifnotfound = list(NA))
    }
    if(!missing(ensembl)) {  # ensembl was supplied
      nms <- ensembl
      eg.list <- mget(ensembl, org.Hs.egENSEMBL2EG, ifnotfound = list(NA))
    }
    .checkLookup(nms, eg.list)
    eg <- sapply(eg.list, function(x) if(length(x) == 1) x else NA)
  }
  if (all(is.na(eg))) {
    out <- rep(NA, length(eg))
  } else { 
        scores <- get(paste("scores.", chip, sep = ""))
        wh <- match(eg, scores$EntrezID)
        out <- rownames(scores)[wh]
        out[is.na(eg)] <- NA
  }
  noProbeset <- nms[is.na(out) & !is.na(eg)]
  if(length(noProbeset) > 0) {
    warning("no probe sets were found for: ", paste(noProbeset, collapse = ", "), call. = FALSE)
  }
  names(out) <- nms
  out
}  

