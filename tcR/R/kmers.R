########## Statistics and analysis of k-mers usage ##########


#' Get kmers from sequences.
#' 
#' @aliases get.kmers
#' 
#' @description
#' Get vector of kmers from the given character vector or data frame.
#' 
#' @param .data Either character vector or a data.frame.
#' @param .head Parameter for head function applied to the given data before kmer generation.
#' @param .k Size of the kmer.
#' @param .clean if T then remove sequences which contain '~' or '*' symbols. Useful for deleting out-of-frame aminoacid sequnces.
#' @param .meat if TRUE than .data must be data.frame with columns CDR3.amino.acid.sequence and Read.count.
#' @param .verbose if T then print progress.
#' @param .left.shift Cut all \code{.left.shift} symbols from the left side for each sequence.
#' @param .right.shift Cut all \code{.right.shift} symbols from the right side for each sequence.
#' 
#' @return Data.frame with 2 columns Kmers and Count / Rank / Proportion relatively to the .value param
#' or a list with such data.frames if .data is a list.
get.kmers <- function (.data, .head = -1, .k = 5, .clean = T, .meat = F, .verbose = T, .left.shift = 0, .right.shift = 0) {  
  if (class(.data) == 'list') {
    ngrams <- lapply(.data, get.kmers, .head = .head, .k = .k, .clean = .clean, .meat = .meat, .left.shift = .left.shift, .right.shift = .right.shift, .verbose = .verbose)
    res <- ngrams[[1]]
    if (.verbose) cat('Merging all data.frames together...\n')
    for (i in 2:length(.data)) {
      cat(i, '/', length(.data), '\n')
      res <- merge(res, ngrams[[i]], by = 'Kmers', all = T)
      names(res) <- c('Kmers', names(.data)[1:i])
    }
    if (.verbose) cat('Done.\n')
    res[is.na(res)] <- 0
    return(res)
  }
  
  .n <- .k  
  
  if (.head == -1) {
    .head <- dim(.data)[1]
    if (is.null(.head)) { .head <- length(.data) }
  }
  
  .data <- head(.data, .head)
  
  read.count <- rep.int(1, .head)
  if (.meat) { read.count <- .data$Read.count }
  if (class(.data) == 'data.frame') { .data <- .data$CDR3.amino.acid.sequence }
  if (.clean) {
    if (.verbose) cat('Cleaning bad sequences...\t')
    .data <- .data[grep('[*, ~]', .data, invert = T)]
    if (.verbose) cat('Sequences after cleaning:', length(.data),'\n')
  }
  
  if (.verbose) cat('Calculating space...\n')
  .data <- substr(.data, 1 + .left.shift, nchar(.data) - .right.shift)
  non.nchar <- nchar(.data) >= .n
  .data <- .data[non.nchar]
  read.count <- read.count[non.nchar]
  space <- sum(nchar(.data) -.n + 1)
  if (.verbose) cat('Number of k-mers:', space,'\n')
  if (space > 0) {
    if (.verbose) {
      cat('Generating k-mers...\n')
      pb <- set.pb(space)
    }
    res <- rep(x='', times=space)
    meat <- rep(1, times = space)
    j <- 1
    for (i in 1:length(.data)) {
      ngrams <- sapply(1:(nchar(.data[i]) - .n + 1), function(j) substr(.data[i], j, j + .n - 1), USE.NAMES=F)
      kmer.meat <- read.count[i]
      for (ngram in ngrams) {
        res[j] <- ngram
        meat[j] <- kmer.meat
        j <- j + 1
        if (.verbose) add.pb(pb)
      }
    }
    if (.verbose) close(pb)
    
    meat <- meat[order(res)]
    res <- sort(res)
    dup <- duplicated(res)
    res.unique <- res[!dup]
    res.dup <- res[dup]
    # res now is a vector for storing number of strings, i.e. res[i] == # of res.unique[i]
    res <- meat[!dup]
    meat <- meat[dup]
  #   res <- rep.int(x=1, times=length(res.unique))
    j <- 1
    if (.verbose) cat('Unique k-mers:', length(res.unique), '\n')
    if (.verbose) cat('Merging k-mers...\n')
    for (i in 1:length(res.unique)) {
      while (res.dup[j] == res.unique[i] && j <= length(res.dup)) {
  #       res[i] <- res[i] + 1
        res[i] <- res[i] + meat[j]
        j <- j + 1
      }
    }
    if (.verbose) cat('Done.\n')
    
    res <- data.frame(Kmers = res.unique, Count = res, stringsAsFactors=F)
    res <- res[order(res$Count, decreasing = T),]
    row.names(res) <- NULL
    res
  
  } else {
    data.frame(Kmers = NA, Count = NA)
  }
}


#' Make and manage the table of the most frequent k-mers.
#' 
#' @aliases kmer.table get.kmer.column
#' 
#' @description
#' \code{kmer.table} - generate table with the most frequent k-mers.
#' 
#' \code{get.kmer.column} - get vector of k-mers from the k-mer table from the function \code{kmer.table}
#' 
#' @usage
#' kmer.table(.data, .heads = c(10, 100, 300, 1000, 3000, 10000, 30000), .nrow = 20,
#'            .clean = T, .meat = F)
#' 
#' get.kmer.column(.kmer.table.list, .head)
#' 
#' @param .data Mitcr data.frame or a list with mitcr data.frames.
#' @param .heads Vector of parameter for the \code{head()} function, supplied sequentialy to the \code{get.kmers()} function. -1 means all rows.
#' @param .nrow How many most frequent k-mers include to the output table.
#' @param .clean Parameter for the \code{get.kmers()} function.
#' @param .meat Parameter for the \code{get.kmers()} function.
#' @param .kmer.table.list Result from the \code{kmer.table} function if \code{.data} supplied as a list.
#' @param .head Which columns with this head return.
#' 
#' @return
#' \code{kmer.table} - if \code{.data} is a data frame, than data frame with columns like "Kmers.X", "Count.X" where X - element from \code{.heads}.
#' If \code{.data} is a list, than list of such data frames.
#' 
#' \code{get.kmer.column} - data frame with first column with kmers and other columns named as a names of data frames, from which \code{.kmer.table.list}
#' was generated.
#' 
#' @examples
#' \dontrun{
#' twb.kmers <- kmer.table(twb, .heads = c(5000, 10000), .meat = T)
#' head(get.kmer.column(twb.kmers, 10000))
#' }
kmer.table <- function (.data, .heads = c(10, 100, 300, 1000, 3000, 10000, 30000), .nrow = 20, .clean=T, .meat = F) {
  if (class(.data) == 'list') {
    return(lapply(.data, kmer.table, .heads = .heads, .nrow = .nrow, .clean = .clean, .meat = .meat))
  }
  res <- do.call(cbind, lapply(.heads, function (h) { head(get.kmers(.data=.data, .head=h, .clean=.clean, .meat = .meat), .nrow) }))
  names(res) <- sapply(1:length(names(res)), function (i) {
    if (.heads[(1+i)%/%2] == -1) { paste0(names(res)[i], '.', 'all', collapse = '') }
    else                         { paste0(names(res)[i], '.', .heads[(1+i)%/%2], collapse = '') }
  })
  res
}

get.kmer.column <- function (.kmer.table.list, .head) {
  name.col <- strsplit(colnames(.kmer.table.list[[1]])[2], '.', fixed=T, useBytes=T)[[1]][1]
  if (.head == -1) { table.names <- paste0(c('Kmers', name.col), '.', 'all') }
  else             { table.names <- paste0(c('Kmers', name.col), '.', .head) }
  kmers <- unique(unlist(lapply(.kmer.table.list, function (x) {x[table.names[1]]})))
  res <- sapply(.kmer.table.list, function (tab) {
    indices <- match(tab[[table.names[1]]], kmers)
    tmp <- rep.int(0, times=length(kmers))
    tmp[indices] <- tab[[table.names[2]]]
    tmp
  })
  df <- data.frame(kmers, res, stringsAsFactors=F)
  names(df) <- c(table.names[1], names(.kmer.table.list))
  df <- df[order(df[, table.names[1]]),]
  row.names(df) <- NULL
  df
}


#' Generate k-mers.
#' 
#' @aliases generate.kmers generate.kmers.prob
#' 
#' @description
#' Generate all k-mers. starting with the given sequence on the given alphabet
#' Generate k-mers with the given k and probabilities of nucleotides next to each other (markov chain).
#' 
#' @usage
#' generate.kmers(.k, .seq = '', .alphabet = c('A', 'C', 'G', 'T'))
#' 
#' generate.kmers.prob(.k, .probs, .count = 1, .alphabet = c('A', 'C', 'G', 'T'),
#'                     .last.nucl = 'X')
#' 
#' @param .k Size of k-mers or either integer or vector with k-s for generate.kmers.prob.
#' @param .seq Prefix of all generated k-mers.
#' @param .probs Matrix with probabilities for generating adjacent symbol with |alphabet| X |alphabet| size. Order of letters is
#' given in the \code{.alphabet} parameter.
#' @param .count Number of kmers to be generated.
#' @param .alphabet Alphabet.
#' @param .last.nucl Adjacent nucleotide from which start generation. If 'X' than choose one of the nucleotides with equal probabilities.
#' 
#' @return Vector of all possible k-mers for \code{generate.kmers}
#' or a vector of generated kmers for \code{generate.kmers.prob}.
generate.kmers <- function (.k, .seq = '', .alphabet = c('A', 'C', 'G', 'T')) {
  kmers.list <- .seq
  for (k in 1:.k) {
    kmers.list <- unlist(lapply(kmers.list, function (kmer) {
      paste0(kmer, .alphabet)
    }))
  }
  kmers.list
}

generate.kmers.prob <- function (.k, .probs, .count = 1, .alphabet = c('A', 'C', 'G', 'T'), .last.nucl = 'X') {  
  if (length(.k) == 1) { .k <- rep.int(.k, .count) }
  if (length(.last.nucl) == 1) { .last.nucl <- rep(.last.nucl, times=.count) }
  
  nucls <- rep('', times = length(.k))
  k.more0 <- .k > 0
  k.more1 <- .k > 1
  nucls[k.more0] <- sapply(.last.nucl[k.more0], function (ln) {
    if (ln == 'X' || ln == '') {
      sample(x=.alphabet, size=1)
    } else {
      sample(x=.alphabet, size=1, prob=.probs[,match(ln, .alphabet)+1], replace = T) 
    }
  })
  nucls[k.more1] <- sapply(which(k.more1), function (i) {  
    nucl <- rep(nucls[i], times = .k[i])
    for (j in 2:.k[i]) {
      nucl[j] <- sample(x=.alphabet, size=1, prob=.probs[,match(nucl[j-1], .alphabet)+1])
    }
    paste0(nucl, collapse='')
  } )
  
  unlist(nucls)
}


#' Profile of sequences of equal length.
#'
#' @aliases kmers.profile
#' 
#' @description
#' Return profile for the given character vector or a data frame with sequences 
#' of equal length or list with them.
#' 
#' @usage
#' kmer.profile(.data, .names = rep('Noname', times=length(.data)), .verbose = F)
#'
#' @param .data Either list with elements of one of the allowed classes or object with one of the class:
#' data.frame with first column with character sequences and second column with number of sequences or character vector.
#' @param .names Names for each sequence / row in the \code{.data}.
#' @param .verbose if T then print processing output.
#' 
#' @return Return data frame with first column "Symbol" with all possible symbols in the given sequences
#' and other columns with names "1", "2", ... for each position with percentage for each symbol.
#' 
#' @seealso \link{vis.logo}
kmer.profile <- function (.data, .names = rep('Noname', times=length(.data)), .verbose = F) {
  .get.nth.letter.stats <- function (.data, .n) {
    res <- summarise(grouped_df(data.frame(Letter = substr(.data[, 1], .n, .n), Count = .data[,2]), list(as.name("Letter"))), Sum.count = sum(Count))
    res$Sum.count <- res$Sum.count / sum(res$Sum.count)
    names(res) <- c("Sequence", as.character(.n))
    res
  }
  
  if (class(.data) == 'list') {
    res <- kmer.profile(.data[[1]], .verbose = .verbose)
    names(res) <- c('Var1', paste0(.names[[1]], c(1,2,3,4,5)))
    for (i in 2:length(.data)) {
      old.names <- names(res)
      res <- merge(res, kmer.profile(.data[[i]]), by = 'Var1', all = T)
      names(res) <- c(old.names, paste0(.names[[i]], c(1,2,3,4,5)))
    }
    return(res)
  }
  if (class(.data) == 'character') {
    .data <- data.frame(Sequence = .data, Read.count = 1, stringsAsFactors=F)
  }
  
  tmp <- .get.nth.letter.stats(.data, 1)
  
  aa.profiles <- tmp
#     aa.profiles <- as.data.frame(prop.table(table(substr(kmers[,1], i, i))), stringsAsFactors=F)
  for (i in 2:nchar(.data[1,1])) {
    if (.verbose) { cat(i, '/', nchar(.data[1,1]), '\n') }
    
    tmp <- .get.nth.letter.stats(.data, i)
    
    aa.profiles <- merge(aa.profiles, tmp, all=T, by = 'Sequence')
#       aa.profiles <- merge(aa.profiles, as.data.frame(prop.table(table(substr(kmers[,1], i, i))), stringsAsFactors=F), all=T, by = 'Var1')
    names(aa.profiles) <- c('Sequence', 1:i)
  }
  aa.profiles <- as.data.frame(aa.profiles)
  names(aa.profiles)[1] <- 'Symbol'
  aa.profiles[is.na(aa.profiles)] <- 0
  aa.profiles <- aa.profiles[order(aa.profiles[,1]),]
  row.names(aa.profiles) <- aa.profiles[,1]
  aa.profiles
}


#' Gibbs Sampler.
#' 
#' @aliases gibbs.sampler
#' 
#' @description
#' Perform the Gibbs Sampler method for finding frequent motifs in the given vector of strings or data.frame.
#' Each string splitted to kmers with the given length of motif.
#' 
#' @param .data Vector of characters or data.frame of characters (1st col) and their numbers (2nd col) if .meat == T.
#' @param .k Motif's length.
#' @param .niter Number of iterations.
#' 
#' @return Vector of possible motifs.
gibbs.sampler <- function (.data, .k = 5, .niter = 500) {
  .n <- .k
  
  .get.best.motif <- function (.seq, .profile) {
    nc <- nchar(.seq)
    res <- rep(x='aaaaa', times=nc - .n + 1)
    for (i in 1:(nc - .n + 1)) {
      res[i] <- substr(.seq, i, i + .n - 1)
    }
    
    max.p <- 0
    max.kmer <- res[1]
    for (km in res) {
      p <- 1
      for (i in 1:.n) {
        tmp <- .profile[match(substr(km, i, i), .profile[,1]),i+1]
        if (!is.na(tmp)) {
          p <- p * tmp
        } else {
          p <- 0
          break
        }
      }
      if (p > max.p) {
        max.p <- p
        max.kmer <- km
      }
    }
    max.kmer
  }
  
  .score <- function (.kmers) {
    sc <- rep.int(0, .n)
    for (i in 1:.n) {
      tab <- table(substr(.kmers, i, i))
      sc[i] <- sum(tab) - max(tab)
    }
    sum(sc)
  }
  
  .data <- .data[nchar(.data) >= .n]
  tmp <- trunc(runif(length(.data), 1, nchar(.data) - .n + 1 + 0.999 ))
  kmers <- substr(.data, tmp, tmp + .n - 1)
  best.kmers <- kmers
  
  cat('Performing gibbs sampling...\n')
  pb <- set.pb(.niter)
  for (i in 1:.niter) {
    unused.motif.index <- sample(1:length(kmers), 1)
    prof <- kmer.profile(kmers[-unused.motif.index])
    kmers[unused.motif.index] <- .get.best.motif(.data[unused.motif.index], prof)
    if (.score(kmers) < .score(best.kmers)) {
      best.kmers <- kmers
    }
    if (i %% 500 == 0) {
      save(best.kmers, file = paste0('gibbs.', i, date(), '.rda'))
    }
    add.pb(pb)
  }
  close(pb)
  
  unique(best.kmers)
}