########## Support functions for managing sequences ##########


# Private function. Split the given string by ', ' to a list and return
# a first element of the list.
.split.get <- function (.str, .alphabet) {
  if (has.class(.str, 'list')) {
    return(lapply(.str, .split.get, .alphabet = .alphabet))
  }
  
  if (has.class(.str, 'data.frame')) {
    res <- .str
    res$V.gene <- sapply(res$V.gene, .split.get, .alphabet = HUMAN_TRBV_MITCR)
    res$D.gene <- sapply(res$D.gene, .split.get, .alphabet = c('TRBD1', 'TRBD2'))
    res$J.gene <- sapply(res$J.gene, .split.get, .alphabet = HUMAN_TRBJ)
    return(res)
  }
  .alphabet2 <- sub(' ', '', .alphabet, fixed = T)
  if (.str %in% .alphabet || .str %in% .alphabet2) {
    .str
  } else {
    #strsplit(.str, ', ', fixed=T, useBytes=T)[[.get]][1]
    strsplit(.str, split="[,][ ]*")[[1]][1]
  }
}


.fix.segments <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, .fix.segments))
  }
  
  .data$V.gene <- sapply(strsplit(.data$V.gene, ',', fixed=T, useBytes=T), paste0, collapse = ', ')
  .data$D.gene <- sapply(strsplit(.data$D.gene, ',', fixed=T, useBytes=T), paste0, collapse = ', ')
  .data$J.gene <- sapply(strsplit(.data$J.gene, ',', fixed=T, useBytes=T), paste0, collapse = ', ')
  .data
}


#' Reverse given character vector by the given n-plets.
#' 
#' @param .seq Sequences.
#' @param .n By which n-plets we should reverse the given strings.
#' 
#' @return Reversed strings.
#' 
#' @examples
#' reverse.string('abcde')  # => "edcba"
#' reverse.string('abcde', 2)  # => "debca"
reverse.string <- function (.seq, .n = 1) {
  lens <- nchar(.seq)
  sapply(1:length(.seq), function (i) {
    tmp <- substring(.seq[i], 1, )
    paste0(c(substring(.seq[i], seq(lens[i] - .n + 1, 1, -.n), 
                       seq(lens[i], .n, -.n)),
             substr(.seq[i], 1, lens[i] %% .n)), 
           collapse = '')
  }, USE.NAMES=F)
}


#' Get all substrings for the given sequence.
#' 
#' @param .seq Sequence for splitting to substrings.
#' @param .min.len Minimal length of output sequences.
#' @param .table if T then return data frame with substrings and positions of their ends in the .seq.
#' 
#' @return Character vector or data frame with columns "Substring", "Start" and "End".
get.all.substrings <- function (.seq, .min.len = 3, .table = T) {
  nc <- nchar(.seq)
  if (.table) {
    seqs <- unlist(sapply(1:(nc - .min.len + 1), function (i) {
      substring(.seq, i, (i + .min.len - 1):nc)
    }))
    inds <- do.call(rbind, sapply(1:(nc - .min.len + 1), function (i) {
      cbind(i, (i + .min.len - 1):nc)
    }))
    data.frame(Substring = seqs, Start = inds[,1], End = inds[,2], stringsAsFactors=F)
  } else {
    unique(unlist(sapply(1:(nc - .min.len + 1), function (i) {
      substring(.seq, i, (i + .min.len - 1):nc)
    })))
  }
}


#' DNA reverse complementing and translation.
#' 
#' @aliases revcomp bunch.translate
#' 
#' @usage
#' revcomp(.seq)
#' 
#' bunch.translate(.seq, .two.way = T)
#' 
#' @description
#' Functions for DNA reverse complementing and translation.
#' 
#' @usage
#' revcomp(.seq)
#' 
#' bunch.translate(.seq, .two.way = T)
#' 
#' @param .seq Vector of nucleotide sequences.
#' @param .two.way if T then translate sequences from both ends (output differes for
#' out-of-frame sequences).
#' 
#' @return Vector of corresponding revese complemented or aminoacid sequences.
revcomp <- function (.seq) {
  rc.table <- c(A = 'T', T = 'A', C = 'G', G = 'C')
  sapply(strsplit(toupper(.seq), '', T, F, T), function (l) paste0(rc.table[l][length(l):1], collapse = ''), USE.NAMES = F)
}

bunch.translate <- function(.seq, .two.way = T) {  
  sapply(toupper(.seq), function (y) {
    ny <- nchar(y)
    ny3 <- ny %/% 3
    tmp <- ''
    if (.two.way) {
      if (ny %% 3 != 0) { tmp <- paste0(rep('N', times = 3), collapse = '') }
      y <- paste0(substr(y, 1, 3*((ny3 %/% 2) +  (ny %% 2))),
                  tmp,
                  substr(y, 3*((ny3 %/% 2) +  (ny3 %% 2)) + (ny %% 3) + 1, ny),
                  collapse = '')
    } else {
      y <- substring(y, seq(1, nchar(y) - 2, 3), seq(3, nchar(y), 3))
    }
    paste0(AA_TABLE[unlist(strsplit(gsub("(...)", "\\1_", y), "_"))],collapse="")
  }, USE.NAMES = F)
}


#' Functions for working with aminoacid sequences.
#' 
#' @aliases codon.variants translated.nucl.sequences reverse.translation
#' 
#' @usage
#' codon.variants(.aaseq, .nucseq = sapply(1:length(.aaseq), 
#'                function (i) paste0(rep('XXX', times = nchar(.aaseq[i])),
#'                collapse = '')))
#' 
#' translated.nucl.sequences(.aaseq, .nucseq = sapply(1:length(.aaseq), 
#'                           function (i) paste0(rep('XXX', times = nchar(.aaseq[i])),
#'                           collapse = '')))
#' 
#' reverse.translation(.aaseq, .nucseq = paste0(rep('XXX', times = nchar(.aaseq)),
#'                    collapse = ''))
#' 
#' @description
#' \code{codon.variants} - get all codon variants for the given nucleotide sequence with known corresponding aminoacid sequence.
#' 
#' \code{translated.nucl.variants} - get number of nucleotide sequences which can be translated to the given aminoacid sequence.
#' 
#' \code{reverse.translation} - get all nucleotide sequences, which can be traslated to the given aminoacid sequence.
#' 
#' @param .aaseq Amino acid sequence.
#' @param .nucseq Nucleotide sequence with 'X' letter at non-fixed positions. Other positions will be fixed.
#' 
#' @return List with all possible variants for every aminoacid in .aaseq, number of sequences or
#' character vector of candidate sequences.
#' 
#' @examples
#' codon.variants('ACT')
#' translated.nucl.sequences(c('ACT', 'CASSLQ'))
#' reverse.translation('T')  # -> "ACA" "ACC" "ACG" "ACT"
#' reverse.translation('T', 'XXT')  # -> "ACT"
#' translated.nucl.sequences('ACT', 'XXXXXXXC')
#' codon.variants('ACT', 'XXXXXXXC')
#' reverse.translation('ACT', 'XXXXXXXC')
codon.variants <- function (.aaseq, .nucseq = sapply(1:length(.aaseq), function (i) paste0(rep('XXX', times = nchar(.aaseq[i])), collapse = ''))) {
  aas <- strsplit(.aaseq, '')
  nuc.nchar <- nchar(.nucseq)
  lapply(1:length(aas), function (i) {
    if (nuc.nchar[i] > 2) {
      triplets <- substring(.nucseq[i], seq(1, nuc.nchar[i] - 2, 3), seq(3, nuc.nchar[i], 3))
      lapply(1:min(length(aas[[i]]), nuc.nchar[i] %/% 3), function (j)
        grep(gsub('X', '[ACGT]', triplets[j]), AA_TABLE_REVERSED[[aas[[i]][j]]], value = T))
    } else {
      list()
    }
  })
}

translated.nucl.sequences <- function (.aaseq, .nucseq = sapply(1:length(.aaseq), function (i) paste0(rep('XXX', times = nchar(.aaseq[i])), collapse = ''))) {
  aas <- strsplit(.aaseq, '')
  nuc.nchar <- nchar(.nucseq)
  sapply(1:length(aas), function (i) {
    if (nuc.nchar[i] > 2) {
      triplets <- substring(.nucseq[i], seq(1, nuc.nchar[i] - 2, 3), seq(3, nuc.nchar[i], 3))
      prod(sapply(1:min(length(aas[[i]]), nuc.nchar[i] %/% 3), function (j)
        length(grep(gsub('X', '[ACGT]', triplets[j]), AA_TABLE_REVERSED[[aas[[i]][j]]], value = T))))
    } else {
      0
    }
  })
}

reverse.translation <- function (.aaseq, .nucseq = paste0(rep('XXX', times = nchar(.aaseq)), collapse = '')) {
  apply(expand.grid(codon.variants(.aaseq, .nucseq)[[1]], stringsAsFactors=F), 1, paste0, collapse = '')
}


#' GC-content of a nucleotide sequences.
#' 
#' @description
#' Compute the GC-content (proportion of G-C nucleotide in a sequence).
#' 
#' @param .nucseq Character vector of nucletoide sequences.
#' @return Numeric vector of \code{length(.nucseq)}.
gc.content <- function (.nucseq) {
  sapply(strsplit(.nucseq, '', fixed=T, useBytes=T), function (l) {
    l <- unlist(l)
    t <- table(c('A', 'C', 'G', 'T', l))
    r <- (t['G'] + t['C'] - 2) / length(l)
    names(r) <- NULL
    r
  }, USE.NAMES=F)
}


#' Find similar sequences.
#' 
#' @aliases find.similar.sequences exact.match hamming.match levenshtein.match
#' 
#' @description
#' Return matrix M with two columns. For each element in row i and column j M[i,j] => distance between pattern(i) and data(j) sequences equal to or less than .max.errors.
#' This function will uppercase .data and remove all strings, which have anything than A-Z letters.
#' 
#' @usage
#' find.similar.sequences(.data, .patterns = c(), .method = c('exact', 'hamm', 'lev'),
#'                        .max.errors = 1, .verbose = T, .clear = F)
#' 
#' exact.match(.data, .patterns = c(), .verbose = T)
#' 
#' hamming.match(.data, .patterns = c(), .max.errors = 1, .verbose = T)
#' 
#' levenshtein.match(.data, .patterns = c(), .max.errors = 1, .verbose = T)
#' 
#' @param .data Vector of strings.
#' @param .patterns Character vector of sequences, which will be used for searching for neighbours.
#' @param .method Which method use: 'exact' for exact matching, 'hamm' for Hamming Distance, 'lev' for Levenshtein distance.
#' @param .max.errors Max Hamming or Levenshtein distance between strings. Doesn't use in 'exact' setting.
#' @param .verbose Should function print progress or not. // DON'T USE IT
#' @param .clear if T then remove all sequences with character "*" or "~".
#' 
#' @return Matrix with two columns [i,j], dist(data(i), data(j)) <= .max.errors.
find.similar.sequences <- function (.data, .patterns = c(), .method = c('exact', 'hamm', 'lev'), .max.errors = 1, .verbose = T, .clear = F) {
  if (.method[1] == 'lev') {
    .clear <- T
  }
  data.old.indices <- grep('[*, ~]', .data, invert = T)
  .data <- toupper(.data[data.old.indices])
  pattern.old.indices <- data.old.indices
  if (length(.patterns) == 0) { 
    .patterns <- .data
  }
  else {
    pattern.old.indices <- 1:length(.patterns)
    if (.clear) {
      pattern.old.indices <- grep('[*, ~]', .patterns, invert = T) 
    }
    .patterns <- toupper(.patterns[pattern.old.indices])
  }
  .fun <- switch(.method[1], exact = .exact_search, hamm = .hamming_search, lev = .levenshtein_search)
  res <- .fun(.data, .patterns, .max.errors, .verbose)
  if (length(res) == 0) {
    matrix(ncol = 2, nrow = 1)
  } else {
    res <- cbind(res[seq(from=1, to=length(res), by=2)], res[seq(from=2, to=length(res), by=2)])
    res <- res[order(res[,1]), ]
    if (is.null(dim(res))) {
      res <- t(as.matrix(res))
    }
    res[,1] <- data.old.indices[res[,1]]
    res[,2] <- pattern.old.indices[res[,2]]
    res
  }
}

.exact_search2 <- function (.data, .patterns, .max.errors, .verbose) {
  ps <- data.frame(P = .patterns, Ind = 1:length(.patterns), stringsAsFactors = F)
  agg <- aggregate(formula = Ind ~ P, data = ps, FUN = function (x) x, simplify = F)
  inds <- lapply(match(.data, agg$P), function (i) 
    if (is.na(i))
      NA
    else
      unlist(agg$Ind[i], use.names = F)
    )
  t(do.call(rbind, sapply(1:length(inds), function (i) 
    if (!is.na(inds[[i]][1]))
      cbind(i, inds[[i]]),
    USE.NAMES = F)))
}

exact.match <- function (.data, .patterns = c(), .verbose = T) {
  find.similar.sequences(.data, .patterns, 'exact', .verbose = .verbose, F)
}

hamming.match <- function (.data, .patterns = c(), .max.errors = 1, .verbose = T) {
  find.similar.sequences(.data, .patterns, 'hamm', .max.errors, .verbose, F)
}

levenshtein.match <- function (.data, .patterns = c(), .max.errors = 1, .verbose = T) {
  find.similar.sequences(.data, .patterns, 'lev', .max.errors, .verbose, T)
}