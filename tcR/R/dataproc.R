if(getRversion() >= "2.15.1")  utils::globalVariables(c("alpha.prob", "beta.prob", "fill_vec", "fill_reads"))


#' Compute the number of deletions in MiTCR data frames.
#' 
#' @aliases get.deletions.alpha get.deletions.beta
#' 
#' @usage
#' get.deletions.alpha(.data, .Vs = segments$TRAV, .Js = segments$TRAJ)
#' 
#' get.deletions.beta(.data, .Vs = segments$TRBV, .Js = segments$TRBJ, .Ds = segments$TRBD)
#' 
#' @description
#' Get deletions for VD, DJ, 5'D and 3'D ends and two columns with 
#' total deletions for VD/DJ and 5'D/3'D deletions for the given mitcr data.frame
#' with 0-indexes columns. Cases, in which deletions cannot be determined, will have -1
#' in their cell.
#' 
#' @param .data Mitcr data.frame.
#' @param .Vs Table of V segments; must have 'V.segment' and 'Nucleotide.sequence' columns.
#' @param .Js Table of J segments; must have 'J.segment' and 'Nucleotide.sequence' columns.
#' @param .Ds Table of D segments; must have 'D.segment' and 'Nucleotide.sequence' columns.
#' 
#' @return Mitcr data.frame with 3 (for alpha chains) or 5 (for beta chains) new columns for deletions.
#' 
#' @details By default, \code{*.table} parameters are taken from the \code{segments} data frame which 
#' can be loaded to your R environment with data(segments). Data for segments has been taken from IMGT.
#' 
#' @examples
#' \dontrun{
#' data(segments)
#' immdata <- get.deletions.beta(.data)
#' immdata.prob <- tcr.prob.df(immdata)
#' }
get.deletions.alpha <- function (.data, .Vs = segments$TRAV, .Js = segments$TRAJ) {
  cat('Preparing data...\n')
  Last.V.nucleotide.position <- .data$Last.V.nucleotide.position
  First.J.nucleotide.position <- .data$First.J.nucleotide.position
  seqs <- .data$CDR3.nucleotide.sequence
  
  cat('Computing V deletions...\n')
  pb <- txtProgressBar(max = length(seqs), style = 3)
  Vs <- sapply(.data$V.gene, function (x) {
    V <- .split.get(x)
    res <- .Vs$Nucleotide.sequence[.Vs$V.alleles == V][1]
    if (is.na(res)) {
      res <- .Vs$Nucleotide.sequence[.Vs$V.alleles == paste(V, '-1', sep = '')][1]
    }
    if (is.na(res)) {
      res <- .Vs$Nucleotide.sequence[.Vs$V.alleles == substr(V, 1, nchar(V) - 2)][1]
    }
    add.pb(pb)
    res
  } )
  V.deletions <- nchar(Vs) - Last.V.nucleotide.position - 1
  close(pb)
  
  cat('Computing J deletions...\n')
  pb <- txtProgressBar(max = length(seqs), style = 3)
  Js <- sapply(.data$J.gene, function (x) {
    J <- .split.get(x)
    res <- .Js$Nucleotide.sequence[.Js$J.alleles == J][1]
    if (is.na(res)) {
      res <- .Js$Nucleotide.sequence[.Js$J.alleles == paste(J, '-1', sep = '')][1]
    }
    if (is.na(res)) {
      res <- .Js$Nucleotide.sequence[.Js$J.alleles == substr(J, 1, nchar(J) - 2)][1]
    }
    add.pb(pb)
    res
  } )
  
  J.deletions <- nchar(Js) - (nchar(seqs) - First.J.nucleotide.position)
  
  cat('Aligning J segments...\n')
  pb <- set.pb(length(Js))
  tmp <- sapply(1:length(Js), function (j) {
    J.seq <- Js[j]
    data.nchar <- nchar(seqs[j])
    .data <- seqs[j]
    J.nchar <- nchar(J.seq)
    add.pb(pb)
  flag <- F
  for (i in 1:(min(data.nchar, J.nchar))) {
    if (substr(J.seq, J.nchar - i + 1, J.nchar - i + 1) != substr(.data, data.nchar - i + 1, data.nchar - i + 1)) {
      flag <- T
      break
    }
  }
  if (flag) i <- i - 1
  substr(J.seq, J.nchar - i + 1, J.nchar)
  })
  close(pb)
  
  J.deletions <- nchar(Js) - nchar(tmp)
  
  close(pb)
  
  .data$V.deletions <- V.deletions
  .data$J.deletions <- J.deletions
  .data$Total.deletions <- V.deletions + J.deletions
  .data$Total.deletions[V.deletions < 0 | J.deletions < 0] <- -1
  .data
}

get.deletions.beta <- function (.data, .Vs = segments$TRBV, .Js = segments$TRBJ, .Ds = segments$TRBD) {
  if (has.class(.data, 'list')) {
    res <- lapply(1:length(.data), function (i) {
      cat('Processing data from', paste0('"', names(.data)[i], '"...', collapse=''), i,'/',length(.data),'\n')
      get.deletions.beta(.data[[i]], .Vs = .Vs, .Js = .Js, .Ds = .Ds)
    } )
    names(res) <- names(.data)
    return(res)
  }
  
  cat('Preparing data...\n')
  pb <- txtProgressBar(max = 5, style = 3)
  Last.V.nucleotide.position <- .data$Last.V.nucleotide.position
  add.pb(pb)
  First.D.nucleotide.position <- .data$First.D.nucleotide.position
  add.pb(pb)
  Last.D.nucleotide.position <- .data$Last.D.nucleotide.position
  add.pb(pb)
  First.J.nucleotide.position <- .data$First.J.nucleotide.position
  add.pb(pb)
  seqs <- .data$CDR3.nucleotide.sequence
  add.pb(pb)
  close(pb)
  
  cat('Computing VD deletions...\n')
  pb <- txtProgressBar(max = length(seqs), style = 3)
  Vs <- sapply(.data$V.gene, function (x) {
    V <- .split.get(x)
    res <- .Vs$Nucleotide.sequence[.Vs$V.segment == V][1]
    if (is.na(res)) {
      res <- .Vs$Nucleotide.sequence[.Vs$V.segment == paste(V, '-1', sep = '')][1]
    }
    if (is.na(res)) {
      res <- .Vs$Nucleotide.sequence[.Vs$V.segment == substr(V, 1, nchar(V) - 2)][1]
    }
    add.pb(pb)
    res
  } )
  VD.deletions <- nchar(Vs) - Last.V.nucleotide.position - 1
  close(pb)
  
  cat('Computing DJ deletions...\n')
  pb <- txtProgressBar(max = length(seqs), style = 3)
  Js <- sapply(.data$J.gene, function (x) {
    J <- .split.get(x)
    res <- .Js$Nucleotide.sequence[.Js$J.segment == J][1]
    if (is.na(res)) {
      res <- .Js$Nucleotide.sequence[.Js$J.segment == paste(J, '-1', sep = '')][1]
    }
    if (is.na(res)) {
      res <- .Js$Nucleotide.sequence[.Js$J.segment == substr(J, 1, nchar(J) - 2)][1]
    }
    add.pb(pb)
    res
  } )
  DJ.deletions <- nchar(Js) - (nchar(seqs) - First.J.nucleotide.position)
  close(pb)
  
  cat("Computing D'5 and D'3 deletions...\n")
  pb <- txtProgressBar(max = length(seqs) * 2, style=3)
  D.gene <- .data$D.gene
  D.straight <- sapply(1:length(D.gene), function (i) {
    add.pb(pb)
    substr(seqs[i],
           First.D.nucleotide.position[i] + 1,
           Last.D.nucleotide.position[i] + 1)
  } )
  D.rev <- revcomp(D.straight)
  D.deletions <- sapply(1:length(D.gene), function (i) {
    D <- .split.get(D.gene[i])
    first.index <- -1
    res <- c(-1, -1)
    if (D != '') {
      D.str <- .Ds$Nucleotide.sequence[.Ds[,1] == D]
      #       D.sub <- substr(seqs[i],
      #                       First.D.nucleotide.position[i] + 1,
      #                       Last.D.nucleotide.position[i] + 1)
      D.sub <- D.straight[i]
      first.index <- regexpr(D.sub, D.str, fixed=T, useBytes=T)[1]
      if (first.index != -1) {
        res[1] <- first.index - 1
        res[2] <- nchar(D.str) - (first.index + nchar(D.sub) - 1)
      } else {
        #         first.index <- regexpr(rev.comp(D.sub), D.str, fixed=T, useBytes=T)[1]
        first.index <- regexpr(D.rev[i], D.str, fixed=T, useBytes=T)[1]
        if (first.index != -1) {
          res[1] <- first.index - 1
          res[2] <- nchar(D.str) - (first.index + nchar(D.sub) - 1)
        }
      }
    }
    add.pb(pb)
    res
  } )  
  close(pb)
  
  .data$VD.deletions <- VD.deletions
  .data$DJ.deletions <- DJ.deletions
  .data$D5.deletions <- D.deletions[1,]
  .data$D3.deletions <- D.deletions[2,]
  .data$Total.deletions <- VD.deletions + DJ.deletions + D.deletions[1,] + D.deletions[2,]
  .data$Total.deletions[VD.deletions < 0 | DJ.deletions < 0 | D.deletions[1,] < 0 | D.deletions[2,] < 0] <- -1
  .data
}


#' Generate random nucleotide TCR sequences.
#' 
#' @description
#' Given the list of probabilities and list of segments (see "Details"), generate a artificial TCR repertoire.
#' 
#' @param .count Number of TCR sequences to generate.
#' @param .chain Either "alpha" or "beta" for alpha and beta chain respectively.
#' @param .segments List of segments (see "Details").
#' @param .P.list List of probabilities (see "Details").
#' 
#' @details
#' For the generation of a artifical TCR repertoire user need to provide two objects: the list with segments and the list with probabilities.
#' List with segments is a list of 5 elements with 5 names: "TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ". Each element is a data frame with following columns 
#' (order is matters!): "V.allelles" with names for V-segments (for TRAV and TRBV; for others is "J.allelles" or "D.allelles"), "CDR3.position" (the function doesn't use it, but you
#' should provide it, fill it with zeros, for example), "Full.nucleotide.sequence" (the function doesn't use it), "Nucleotide.sequence" (function uses it for getting nucleotide
#' sequences of segments) and "Nucleotide.sequence.P" (the function doesn't use it).
#' 
#' List with probabilities is quite complicated, so just call \code{data(beta.prob)} for beta chain probabilities (alpha chain probabilities will be added soon).
#' 
#' @return Mitcr data.frame with generated sequences.
#' 
#' @seealso \link{genesegments} \link{beta.prob}
#' 
#' @examples
#' \dontrun{
#' # Load list of segments provided along with tcR.
#' data(genesegments) 
#' # Load list of probabilities provided along with tcR.
#' data(beta.prob)
#' # Generate repertoire of beta chian with 10000 sequences.
#' artif.rep <- generate.tcR(10000, 'beta')
#' View(artif.rep)
#' }
generate.tcr <- function (.count = 1, .chain = c('beta', 'alpha'), .segments,
                          .P.list = if (.chain[1] == 'alpha') alpha.prob else beta.prob) {
  #   .unique = -1
  
  cat('Preparing necessary data...\n')
  ###############
  # ALPHA CHAIN #
  ###############
  if (.chain[1] == 'alpha') {
    pb <- set.pb(3)
    .V.J.prob <- .P.list$P.V.J
    V.J.sample <- sample2D(.V.J.prob, .count)
    V.sample <- V.J.sample[,1]
    add.pb(pb)
    J.sample <- V.J.sample[,2]
    add.pb(pb)
    .del.V.prob <- .P.list$P.del.V
#     names(.del.V.prob)[-1] <- sapply(names(.del.V.prob)[-1], function (x) { sub('.', '-', x, fixed=T) })
    V.dels <- sapply(V.sample, function (x) { sample(x=0:(nrow(.del.V.prob) - 1), size=1, prob=.del.V.prob[,x]) } )
    add.pb(pb)
    .del.J.prob <- .P.list$P.del.J
#     names(.del.J.prob)[-1] <- sapply(names(.del.J.prob)[-1], function (x) { sub('.', '-', x, fixed=T) })
    J.dels <- sapply(J.sample, function (x) { sample(x=0:(nrow(.del.J.prob) - 1), size=1, prob=.del.J.prob[,x]) } )
    add.pb(pb)
    
    cat('Generating segments and insertions...\n')
    pb <- set.pb(4)
    .Vs <- .segments$TRAV
    V <- .Vs$Nucleotide.sequence[match(V.sample, .Vs[,1])]
    add.pb(pb)
    .Js <- .segments$TRAJ
    J <- .Js$Nucleotide.sequence[match(J.sample, .Js[,1])]
    add.pb(pb)

    .ins.len.prob <- .P.list$P.ins.len
    VJ.ins.lens <- sample(x=0:(nrow(.ins.len.prob) - 1), size=.count, replace=T, prob=.ins.len.prob[,2])
    add.pb(pb)
    
    .ins.nucl.prob <- .P.list$P.ins.nucl
    VJ.ins.nucls <- generate.kmers.prob(VJ.ins.lens, .probs=.ins.nucl.prob[,c(1, 2:5)], .last.nucl = substr(V, nchar(V) - V.dels, nchar(V) - V.dels))
    add.pb(pb)
    close(pb)
    
    cat('Finalising sequences...\n')
    pb <- set.pb(.count)
    Vns <- nchar(V)
    Jns <- nchar(J)
    seqs <- sapply(1:.count, function (i) { 
      add.pb(pb)
      paste0(substr(V[i], 1, Vns[i] - V.dels[i]),
             VJ.ins.nucls[i],
             substr(J[i], J.dels[i] + 1, Jns[i]),
             collapse = '')
      
    } )
    close(pb)
    
    cat('Generating data frame...\n')
    res<-data.frame(Read.count = 1,
                    CDR3.nucleotide.sequence = seqs,
                    CDR3.amino.acid.sequence = bunch.translate(seqs),
                    V.gene = V.sample,
                    J.gene = J.sample,
                    D.gene = '',
                    VD.insertions = -1,
                    DJ.insertions = -1,
                    Total.insertions = VJ.ins.lens,
                    VJ.inserted.sequence = VJ.ins.nucls,
                    VD.deletions = V.dels,
                    DJ.deletions = J.dels,
                    Total.deletions = V.dels + J.dels,
                    stringsAsFactors=F)
#     if (.inframe){
#       res<-cbind(rep(1,nrow(res)),res)
#       names(res)[1]<-"Read.count"
#       res<-get.outframes((res),T)
#     }
  }
  ##############
  # BETA CHAIN #
  ##############
  else {
    pb <- set.pb(7)
    .V.prob <- .P.list$P.V
    V.sample <- sample(x=row.names(.V.prob), size=.count, replace=T, prob=.V.prob)
    add.pb(pb)
    .D.J.prob <-.P.list$P.J.D
    JD.sample <- sample2D(.D.J.prob, .count=.count)
    add.pb(pb)
    .del.V.prob <- .P.list$P.del.V
#     names(.del.V.prob)[-1] <- sapply(names(.del.V.prob)[-1], function (x) { sub('.', '-', x, fixed=T) })
    VD.dels <- sapply(V.sample, function (x) { sample(x=0:(nrow(.del.V.prob) - 1), size=1, prob=.del.V.prob[,x]) } )
    add.pb(pb)
    .del.J.prob <- .P.list$P.del.J
#     names(.del.J.prob)[-1] <- sapply(names(.del.J.prob)[-1], function (x) { sub('.', '-', x, fixed=T) })
    DJ.dels <- apply(JD.sample, 1, function (x) { sample(x=0:(nrow(.del.J.prob) - 1), size=1, prob=.del.J.prob[,x[1]]) } )
    add.pb(pb)
    
    D1s <- JD.sample[,2] == 'TRBD1'
    D5.D3.dels <- matrix(rep('', times=length(JD.sample)), ncol = 2)
    add.pb(pb)
    .del.D1.prob <- .P.list$P.del.D1
    row.names(.del.D1.prob) <- 0:(nrow(.del.D1.prob) - 1)
    D5.D3.dels[D1s,] <- sample2D(.del.D1.prob, .count=length(which(D1s)))
    add.pb(pb)
    .del.D2.prob <- .P.list$P.del.D2
    row.names(.del.D2.prob) <- 0:(nrow(.del.D2.prob) - 1)
    D5.D3.dels[!D1s,] <- sample2D(.del.D2.prob, .count=dim(JD.sample)[1] - length(which(D1s)))
    add.pb(pb)
    
    cat('Generating segments and insertions...\n')
    pb <- set.pb(7)
    .Vs <- .segments$TRBV
    V <- .Vs$Nucleotide.sequence[match(V.sample, .Vs[,1])]
    add.pb(pb)
    .Js <- .segments$TRBJ
    J <- .Js$Nucleotide.sequence[match(JD.sample[,1], .Js[,1])]
    add.pb(pb)
    .Ds <- .segments$TRBD
    D <- .Ds$Nucleotide.sequence[match(JD.sample[,2], .Ds[,1])]
    add.pb(pb)

    .ins.len.prob <- .P.list$P.ins.len
    VD.ins.lens <- sample(x=0:(nrow(.ins.len.prob) - 1), size=.count, replace=T, prob=.ins.len.prob[,2])
    add.pb(pb)
    DJ.ins.lens <- sample(x=0:(nrow(.ins.len.prob) - 1), size=.count, replace=T, prob=.ins.len.prob[,3])
    add.pb(pb)

    # Few secs here
    .ins.nucl.prob <- .P.list$P.ins.nucl
    VD.ins.nucls <- generate.kmers.prob(VD.ins.lens, .probs=.ins.nucl.prob[,c(1, 2:5)], .last.nucl = substr(V, nchar(V) - VD.dels, nchar(V) - VD.dels))
    add.pb(pb)
    # Few secs here
    DJ.ins.nucls <- generate.kmers.prob(DJ.ins.lens, .probs=.ins.nucl.prob[,c(1, 6:9)], .last.nucl = substr(D, nchar(D) - as.integer(D5.D3.dels[,2]), nchar(D) - as.integer(D5.D3.dels[,2])))
    add.pb(pb)
    close(pb)
    close(pb)
    
    cat('Finalising sequences...\n')
    pb <- set.pb(.count)
    Vns <- nchar(V)
    Dns <- nchar(D)
    Jns <- nchar(J)
    seqs <- sapply(1:.count, function (i) { 
      add.pb(pb)
      paste0(substr(V[i], 1, Vns[i] - VD.dels[i]),
             VD.ins.nucls[i],
             substr(D[i], as.integer(D5.D3.dels[i,1]) + 1, Dns[i] - as.integer(D5.D3.dels[i,2])),
             DJ.ins.nucls[i],
             substr(J[i], DJ.dels[i] + 1, Jns[i]),
             collapse = '')
      
    } )
    seqs.sep <- sapply(1:.count, function (i) { 
      add.pb(pb)
      paste(substr(V[i], 1, Vns[i] - VD.dels[i]),
             VD.ins.nucls[i],
             substr(D[i], as.integer(D5.D3.dels[i,1]) + 1, Dns[i] - as.integer(D5.D3.dels[i,2])),
             DJ.ins.nucls[i],
             substr(J[i], DJ.dels[i] + 1, Jns[i]),
             sep = ':')
      
    } )
    close(pb)
    
    cat('Generating data frame...\n')
    res<-data.frame(Read.count = 1,
                    CDR3.nucleotide.sequence = seqs,
                    CDR3.nucleotide.sequence.sep = seqs.sep,
                    CDR3.amino.acid.sequence = bunch.translate(seqs),
                    V.gene = V.sample,
                    J.gene = JD.sample[,1],
                    D.gene = JD.sample[,2],
                    VD.insertions = VD.ins.lens,
                    DJ.insertions = DJ.ins.lens,
                    Total.insertions = VD.ins.lens + DJ.ins.lens,
                    VJ.inserted.sequence = VD.ins.nucls,
                    DJ.inserted.sequence = DJ.ins.nucls,
                    VD.deletions = VD.dels,
                    DJ.deletions = DJ.dels,
                    Total.deletions = VD.dels + DJ.dels,
                    stringsAsFactors=F)
#     if (.inframe){
#       res<-cbind(rep(1,nrow(res)),res)
#       names(res)[1]<-"Read.count"
#       res<-get.outframes((res),T)
#     }
  }
  
  res$Inframe <- nchar(res$CDR3.nucleotide.sequence) %% 3 == 0
  cat('Generated', sum(res$Inframe), 'in-frame sequences and', nrow(res) - sum(res$Inframe), 'out-of-frame sequences.')
  res
}


#' Proportions of specifyed subsets of clones.
#' 
#' @aliases tailbound.proportion clonal.proportion top.proportion
#' 
#' @usage
#' tailbound.proportion(.data, .bound = 2, .col = 'Read.count')
#' 
#' top.proportion(.data, .head = 10, .col = 'Read.count')
#' 
#' clonal.proportion(.data, .perc = 10, .col = 'Read.count')
#' 
#' @description
#' Get a specifyed subset of the given data and compute which proportion in counts
#' it has comparing to the overall count.
#' 
#' \code{tailbound.proportion} - subset by the count;
#' 
#' \code{top.proportion} - subset by rank (top N clones);
#' 
#' \code{clonal.proportion} - subset by a summary percentage (top N clones which in sum has the given percentage).
#' 
#' @param .data Data frame or a list with data frames.
#' @param .bound Subset the \code{.data} by \code{.col} <= \code{.bound}.
#' @param .head How many top values to choose - parameter to the \code{.head} function.
#' @param .perc Percentage (0 - 100).
#' @param .col Column's name with counts of sequences.
#' 
#' @return For \code{tailbound.proportion} - numeric vector of percentage.
#' 
#' For \code{top.proportion} - numeric vector of percentage for top clones.
#' For \code{clonal.proportion} - vector or matrix with values for number of clones, occupied percentage and proportion of the
#' chosen clones to the overall count of clones.
#' @seealso \link{vis.top.proportions}
#' 
#' @examples
#' \dontrun{
#'                                 # How many clones fill up approximately
#' clonal.proportion(immdata, 25)  # the 25% of the sum of values in 'Read.count'?
#' 
#'                                 # What proportion of the top-10 clones' reads
#" top.proportion(immdata, 10)    # to the overall number of reads?
#' vis.top.proportions(immdata)  # Plot this proportions.
#' 
#'                                 # What proportion of sequences which
#'                                 # has 'Read.count' <= 100 to the
#' tailbound.proportion(immdata, 100)  # overall number of reads?
#' }
tailbound.proportion <- function (.data, .bound = 2, .col = 'Read.count') {
  if (has.class(.data, 'list')) {
    return(sapply(.data, tailbound.proportion, .bound = .bound, .col = .col))
  }
  sum(.data[, .col] <= .bound) / nrow(.data)
}

top.proportion <- function (.data, .head = 10, .col = 'Read.count') {
  if (has.class(.data, 'list')) {
    return(sapply(.data, top.proportion, .head = .head, .col = .col))
  }
  sum(head(.data[, .col], .head)) / sum(.data[, .col])
}

clonal.proportion <- function (.data, .perc = 10, .col = 'Read.count') {
  if (has.class(.data, 'list')) {
    return(t(sapply(.data, clonal.proportion, .perc = .perc, .col = .col)))
  }
  prop <- 0
  n <- 0
  col <- .data[, .col]
  col.sum <- sum(col)
  while (prop < col.sum * (.perc / 100)) {
    n <- n + 1
    prop <- prop + col[n]
  }
  c(Clones = n, Percentage = 100 * signif((prop / col.sum), 3), Clonal.count.prop = n / nrow(.data))
}


#' Rearrange columns with counts for clonesets.
#' 
#' @description
#' Replace Read.count with Umi.count, recompute Percentage and sort data.
#' 
#' @param .data Data frame with columns "Umi.count" and "Read.count".
#' 
#' @return Data frame with new "Read.count" and "Percentage" columns.
barcodes.to.reads <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, barcodes.to.reads))
  }
  .data$Read.count <- .data$Umi.count
  .data$Percentage <- .data$Read.count / sum(.data$Read.count)
  .data[order(.data$Percentage, decreasing = T),]
}


#' Resample data frame using values from the column with number of clonesets.
#' 
#' @aliases resample downsample
#' 
#' @description
#' Resample data frame using values from the column with number of clonesets. Number of clonestes (i.e., rows of a MiTCR data frame)
#' are reads (usually the "Read.count" column) or UMIs (i.e., barcodes, usually the "Umi.count" column).
#' 
#' @param .data Data frame with the column \code{.col} or list of such data frames.
#' @param .n Number of values / reads / UMIs to choose.
#' @param .col Which column choose to represent quanitites of clonotypes. See "Details".
#' 
#' @return Data frame with \code{sum(.data[, .col]) == .n}.
#' 
#' @details
#' \code{resample}. Using multinomial distribution, compute the number of occurences for each cloneset, than remove zero-number clonotypes and
#' return resulting data frame. Probabilities for \code{rmultinom} for each cloneset is a percentage of this cloneset in
#' the \code{.col} column. It's a some sort of simulation of how clonotypes are chosen from the organisms. For now it's not working
#' very well, so use \code{downsample} instead.
#' 
#' \code{downsample}. Choose \code{.n} clones (not clonotypes!) from the input repertoires without any probabilistic simulation, but
#' exactly computing each choosed clones. Its output is same as for \code{resample} (repertoires), but is more consistent and
#' biologically pleasant.
#' 
#' @seealso \link{rmultinom}
#' 
#' @examples
#' \dontrun{
#' # Get 100K reads (not clones!).
#' immdata.1.100k <- resample(immdata[[1]], 100000, .col = "read.count")
#' }
resample <- function (.data, .n = -1, .col = 'read.count') {
  if (has.class(.data, 'list')) {
    if (length(.n) != length(.data)) {
      .n <- c(.n, rep.int(-1, length(.data) - length(.n)))
    }
    return(lapply(.data, resample, .n = .n, .col = .col))
  }
  .col <- .column.choice(.col[1])
  if (.n == -1) {
    .n <- sum(.data[, .col])
  }
  new.bc <- rmultinom(1, .n, .data[, .col] / sum(.data[, .col]))
  .data[, .col] <- new.bc
  non.zeros <- new.bc != 0
  .data <- .data[non.zeros, ]
  perc.col <- paste0(strsplit(.col, ".", T)[[1]][1], ".proportion")
  .data[[perc.col]] <- new.bc[non.zeros] / sum(new.bc)
  .data[order(.data[[perc.col]], decreasing = T),]
}

downsample <- function (.data, .n, .col = c("read.count", "umi.count")) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, downsample, .n = .n, .col = .col))
  }
  
  col_current <- .column.choice(.col[1])
  read_vec <- .data[, col_current]
  read_indices <- rep(0, sum(read_vec))
  cppFunction(
    '
  NumericVector fill_vec(NumericVector read_vec, NumericVector read_indices) 
  {
    int dummy = 0;
    for (int i = 0; i < read_vec.size(); i++)
    {
      for (int j = dummy; j < (read_vec[i] + dummy); j++)
  	{
	    read_indices[j] = i;
  	}
  	dummy = dummy + read_vec[i];
    }
    return read_indices;
  }
  '
  )
  read_indices <- fill_vec(read_vec, read_indices)
  new_counts <- sample(read_indices, .n)
  new_reads <- rep(0, length(read_vec))
  cppFunction(
    '
  NumericVector fill_reads(NumericVector new_reads, NumericVector new_counts) 
  {
    for (int i = 0; i < new_counts.size(); i++)
    {
	    new_reads[new_counts[i]] = new_reads[new_counts[i]] + 1;
	  }
    return new_reads;
  }
  '
  )
  .data[, col_current] <- fill_reads(new_reads, new_counts)
  
  subset(.data, .data[, col_current] > 0)
}


#' Set new columns "Rank" and "Index".
#' 
#' @aliases set.rank set.index
#' 
#' @description
#' Set new columns "Rank" and "Index":
#' 
#' set.rank <==> .data$Rank = rank(.data[, .col], ties.method = 'average')
#' 
#' set.index <==> .data$Index = 1:nrow(.data) in a sorted data frame by \code{.col}
#' 
#' @param .data Data frame or list with data frames.
#' @param .col Character vector with name of the column to use for ranking or indexing.
#' 
#' @return Data frame with new column "Rank" or "Index" or list with such data frames.
set.rank <- function (.data, .col = "Read.count") {
  if (has.class(.data, 'list')) {
    lapply(.data, set.rank, .col = .col)
  } else {
    y <- 1 / .data[[.col]]
    .data$Rank <- rank(y, ties.method = 'average')
    .data
  }
}

set.index <- function (.data, .col = "Read.count") {
  if (has.class(.data, 'list')) {
    return(lapply(.data, set.index))
  }
  .data <- .data[order(.data[, .col], decreasing = T), ]
  .data$Index <- 1:nrow(.data)
  .data
}