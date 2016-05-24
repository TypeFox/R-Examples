########## Repertoire statistics ##########


#' In-frame / out-of-frame sequences filter.
#' 
#' @aliases get.inframes get.outframes count.inframes count.outframes get.frames count.frames clonotypescount
#' 
#' @description
#' Return the given data frame with in-frame or out-of-frame sequences only. Nucleotide sequences in a column "CDR3.nucleotide.sequence" are checked if
#' they length are divisible by 3 (len mod 3 == 0 => in-frame, else out-of-frame)
#' 
#' @usage
#' get.inframes(.data, .head = 0, .coding = T)
#' 
#' get.outframes(.data, .head = 0)
#' 
#' count.inframes(.data, .head = 0, .coding = T)
#' 
#' count.outframes(.data, .head = 0)
#' 
#' get.frames(.data, .frame = c('in', 'out', 'all'), .head = 0, .coding = T)
#' 
#' count.frames(.data, .frame = c('in', 'out', 'all'), .head = 0, .coding = T)
#' 
#' @param .data MiTCR data.frame or a list with mitcr data.frames.
#' @param .frame Which *-frames to choose.
#' @param .head Parameter to the head() function. Supply 0 to get all elements. \code{head} applied before subsetting, i.e.
#' if .head == 500, you will get in-frames from the top 500 clonotypes.
#' @param .coding if T then return only coding sequences, i.e. without stop-codon.
#' 
#' @return Filtered data.frame or a list with such data.frames.
get.inframes <- function (.data, .head = 0, .coding = T) { 
  if (class(.data) == 'list') { return(lapply(.data, get.inframes, .head = .head, .coding = .coding)) }
  .data <- head(.data, if (.head == 0) {nrow(.data)} else {.head})
  if (.coding) {
    d <- subset(.data, nchar(.data$CDR3.nucleotide.sequence) %% 3 == 0)
    d[grep('[*, ~]', d$CDR3.amino.acid.sequence, invert = T), ]
  } else {
    subset(.data, nchar(.data$CDR3.nucleotide.sequence) %% 3 == 0)
  }
}

get.outframes <- function (.data, .head = 0) {
  if (class(.data) == 'list') { return(lapply(.data, get.outframes, .head = .head)) }
  .data <- head(.data, if (.head == 0) {nrow(.data)} else {.head})
  subset(.data, nchar(.data$CDR3.nucleotide.sequence) %% 3 != 0)
}

count.inframes <- function (.data, .head = 0, .coding = T) {
  if (class(.data) == 'list') { sapply(get.inframes(.data, .head, .coding), nrow) }
  else { nrow(get.inframes(.data, .head, .coding)) }
}

count.outframes <- function (.data, .head = 0) {
  if (class(.data) == 'list') { sapply(get.outframes(.data, .head), nrow) }
  else { nrow(get.outframes(.data, .head)) }
}

get.frames <- function (.data, .frame = c('in', 'out', 'all'), .head = 0, .coding = T) {
  if (.frame[1] == 'in') { get.inframes(.data, .head, .coding) }
  else if (.frame[1] == 'out') { get.outframes(.data, .head) }
  else { head(.data, if (.head == 0) {nrow(.data)} else {.head}) }
}

count.frames <- function (.data, .frame = c('in', 'out', 'all'), .head = 0, .coding = T) {
  if (.frame[1] == 'in') { count.inframes(.data, .head, .coding) }
  else if (.frame[1] == 'out') { count.outframes(.data, .head) }
  else { nrow(head(.data, if (.head == 0) {nrow(.data)} else {.head})) }
}

clonotypescount <- function(.data, .head = 0) {
  length(unique(as.character(head(.data, if (.head == 0) {nrow(.data)} else {.head})$CDR3.amino.acid.sequence)))
}


#' MiTCR data frame basic statistics.
#' 
#' @aliases cloneset.stats repseq.stats
#' 
#' @usage
#' cloneset.stats(.data, .head = 0)
#' 
#' repseq.stats(.data, .head = 0)
#' 
#' @description
#' Compute basic statistics of TCR repertoires: number of clones, number of clonotypes, 
#' number of in-frame and out-of-frame sequences, summary of "Read.count", "Umi.count" and other.
#' 
#' @param .data tcR data frames or a list with tcR data frames.
#' @param .head How many top clones use to comput summary.
#' 
#' @return if \code{.data} is a data frame, than numeric vector with statistics. If \code{.data} is 
#' a list with data frames, than matrix with statistics for each data frame.
#' 
#' @examples
#' \dontrun{
#' # Compute basic statistics of list with data frames.
#' cloneset.stats(immdata)
#' repseq.stats(immdata)
#' }
cloneset.stats <- function (.data, .head = 0) {
  if (has.class(.data, 'list')) {
    res <- t(do.call(cbind, lapply(.data, cloneset.stats, .head = .head)))
    row.names(res) <- names(.data)
    return(res)
  }
  .head <- if (.head == 0) {nrow(.data)} else {.head}  
  res <- c(as.numeric(.head),
           clonotypescount(.data, .head),
           clonotypescount(.data, .head) / .head,
           count.inframes(.data, .head),
           count.inframes(.data, .head) / .head,
           count.outframes(.data, .head),
           count.outframes(.data, .head) / .head)
  names(res) <- c('#Nucleotide clones','#Aminoacid clonotypes', '%Aminoacid clonotypes', '#In-frames', '%In-frames', '#Out-of-frames', '%Out-of-frames')
  
  .data <- head(.data, .head)
  res2 <- c(Sum = sum(.data$Read.count), summary(.data$Read.count))
  names(res2) <- sub('.', '', names(res2), fixed = T)
  names(res2) <- paste0(names(res2), '.reads')
  if (!is.na(.data$Umi.count)[1]) {
    res3 <- c(Sum = sum(.data$Umi.count), summary(.data$Umi.count))
    names(res3) <- sub('.', '', names(res3), fixed = T)
    names(res3) <- paste0(names(res3), '.UMIs')
    res2 <- c(res2, res3)
  }
  c(res, res2)
}

repseq.stats = function(.data, .head=0) {
  if (has.class(.data, "list")) {
    res=do.call(cbind, lapply(.data, repseq.stats, .head=.head))
    dimnames(res)[[2]]=names(.data)
    return(t(res))
  }else{
    .umi <- !is.na(.data$Umi.count[1])
    .head= if (.head==0){nrow(.data)} else {.head}
    .data=head(.data, .head)
    if (!.umi) {
      res=c(nrow(.data), sum(.data$'Read.count'), round(sum(.data$'Read.count')/nrow(.data), digits = 2))
      names(res)=c('Clones', "Sum.reads", "Reads.per.clone")
      return(res)
    }else{
      res=c(nrow(.data), sum(.data$'Read.count'), sum(.data$'Umi.count'), round(sum(.data$'Read.count') / sum(.data$'Umi.count'), digits = 2), round(sum(.data$'Umi.count') / nrow(.data), digits = 2))
      names(res)=c('Clones', "Sum.reads", "Sum.UMIs", "Reads.per.UMI", "UMI.per.clone")
      return(res)
    }
  }
}


#' Columns statistics.
#' 
#' @aliases column.summary insertion.stats
#' 
#' @usage
#' column.summary(.data, .factor.col, .target.col, .alphabet = NA, .remove.neg = T)
#' 
#' insertion.stats(.data)
#' 
#' @description
#' \code{column.summary} - general function for computing summary statistics (using the \code{summary} function) for columns of the given mitcr data.frame:
#' divide \code{.factor.column} by factors from \code{.alphabet} and compute statistics
#' of correspondingly divided \code{.target.column}.
#' 
#' \code{insertion.stats} - compute statistics of insertions for the given mitcr data.frame.
#' 
#' @param .data Data frame with columns \code{.factor.col} and \code{target.col}
#' @param .factor.col Columns with factors by which the data will be divided to subsets.
#' @param .target.col Column with numeric values for computing summaries after dividing the data to subsets.
#' @param .alphabet Character vector of factor levels. If NA than use all possible elements frim the \code{.factor.col} column.
#' @param .remove.neg Remove all elements which >-1 from the \code{.target.col} column.
#' 
#' @seealso \link{summary}
#' 
#' @return Data.frame with first column with levels of factors and 5 columns with output from the \code{summary} function.
#' 
#' @examples
#' \dontrun{
#' # Compute summary statistics of VD insertions
#' # for each V-segment using all V-segments in the given data frame.
#' column.summary(immdata[[1]], 'V.gene', 'Total.insertions')
#' # Compute summary statistics of VD insertions for each V-segment using only V-segments
#' # from the HUMAN_TRBV_MITCR
#' column.summary(immdata[[1]], 'V.gene', 'Total.insertions', HUMAN_TRBV_MITCR)
#' }
column.summary <- function (.data, .factor.col, .target.col, .alphabet = NA, .remove.neg = T) {
  if (length(.alphabet) != 0 && !is.na(.alphabet[1])) {
    .data[!(.data[, .factor.col] %in% .alphabet), .factor.col] <- 'Other'
  }

  if (.remove.neg) {
    .data <- .data[.data[, .target.col] > -1, ]
  }
  
  res <- do.call(rbind, tapply(.data[, .target.col], .data[, .factor.col], function (x) c(summary(x))))
  res <- data.frame(A = row.names(res), res, stringsAsFactors=F)
  names(res)[1] <- .factor.col
  row.names(res) <- NULL
  res
}

insertion.stats <- function (.data) {
  if (class(.data) == 'list') {
    return(lapply(.data, insertion.stats))
  }
  vd <- column.summary(.data, 'V.gene', 'VD.insertions', HUMAN_TRBV_MITCR)
  names(vd)[-1] <- paste0('VD.', names(vd)[-1])
  dj <- column.summary(.data, 'V.gene', 'DJ.insertions', HUMAN_TRBV_MITCR)
  names(dj)[-1] <- paste0('DJ.', names(dj)[-1])
  tot <- column.summary(.data, 'V.gene', 'Total.insertions', HUMAN_TRBV_MITCR)
  names(tot)[-1] <- paste0('Total.', names(tot)[-1])
  res <- merge(vd, dj, by = 'V.gene', all = T)
  res <- merge(res, tot, by = 'V.gene', all = T)
  res
}


#' Find target clonotypes and get columns' value corresponded to that clonotypes.
#' 
#' @description
#' Find the given target clonotypes in the given list of data.frames and get corresponding values of desired columns.
#' 
#' @param .data List with mitcr data.frames or a mitcr data.frame.
#' @param .targets Target sequences or elements to search. Either character vector or a matrix / data frame (not a data table!) with two columns: first for sequences, second for V-segments.
#' @param .method Method, which will be used to find clonotypes:
#' 
#' - "exact" performs exact matching of targets;
#' 
#' - "hamm" finds targets and close sequences using hamming distance <= 1;
#' 
#' - "lev" finds targets and close sequences using levenshtein distance <= 1.
#' 
#' @param .col.name Character vector with column names which values should be returned.
#' @param .target.col Character vector specifying name of columns in which function will search for a targets.
#' Only first column's name will be used for matching by different method, others will match exactly.
#' \code{.targets} should be a two-column character matrix or data frame with second column for V-segments.
#' @param .verbose if T then print messages about the search process.
#' 
#' @return Data.frame.
#' 
#' @examples
#' \dontrun{
#' # Get ranks of all given sequences in a list of data frames.
#' immdata <- set.rank(immdata)
#' find.clonotypes(.data = immdata, .targets = head(immdata[[1]]$CDR3.amino.acid.sequence),
#'                 .method = 'exact', .col.name = "Rank", .target.col = "CDR3.amino.acid.sequence")
#' # Find close by levenhstein distance clonotypes with similar V-segments and return
#' # their values in columns 'Read.count' and 'Total.insertions'.
#' find.clonotypes(.data = twb, .targets = twb[[1]][, c('CDR3.amino.acid.sequence', 'V.gene')],
#'                 .col.name = c('Read.count', 'Total.insertions'), .method = 'lev',
#'                 .target.col = c('CDR3.amino.acid.sequence', 'V.gene'))
#' }
find.clonotypes <- function (.data, .targets, .method = c('exact', 'hamm', 'lev'), .col.name = 'Read.count', .target.col = 'CDR3.amino.acid.sequence', .verbose = T) {  
  if (is.character(.targets) && length(.target.col) != 1) {
    cat("Target columns doesn't match the given .targets!\n")
    return()
  }
  
  if (!is.character(.targets) && ncol(.targets) != length(.target.col)) {
    cat("Target columns doesn't match the given .targets!\n")
    return()
  }
  
  if (has.class(.data, 'data.frame')) {
    .data <- list(Sample = .data)
  }
  
  if (.method[1] == 'lev') {
    .data <- lapply(.data, function (.data) .data[grep('[*, ~]', .data[, .target.col[1]], invert = T),])
  }
  
  .verbose.msg('Searching for targets...\n', .verbose)
  if (.verbose) pb <- set.pb(length(.data))
  dt.list <- list()
  for (i in 1:length(.data)) {
    if (.verbose) add.pb(pb)
    if (length(.target.col) ==  1) {
      inds <- intersectIndices(.targets, .data[[i]][, .target.col], .method)
    } else {
      colnames(.targets) <- .target.col
      inds <- intersectIndices(.targets, .data[[i]][, .target.col], .method, .target.col)
    }

    inds <- inds[!duplicated(inds[,2]),]
    
    if(!is.matrix(inds)) {
      inds <- rbind(inds)
    }
    
    if (dim(inds)[1] == 0) {
      inds <- matrix(c(NA, NA), 1, 2)
    }

    res <- list()
    
    if (!is.na(inds[1,1])) {      
      for (col.name in .col.name) {
        res[[col.name]] <- .data[[i]][inds[,2], col.name]
      }
      dt.list[[names(.data)[i]]] <- as.data.table(data.frame(.data[[i]][inds[,2], .target.col], res, stringsAsFactors = F), F)
      setnames(dt.list[[names(.data)[i]]], c(.target.col, .col.name))
      setkeyv(dt.list[[names(.data)[i]]], .target.col)
      tc <- dt.list[[names(.data)[i]]][[.target.col[1]]]
      dups <- duplicated(tc)
      dups.i <- rep.int(0, length(dups))
      k <- 1
      for (j in 1:length(dups)) {
        if (dups[j]) {
          dups.i[j] <- k
          k <- k + 1
        } else {
          k <- 1
        }
      }
      dt.list[[names(.data)[i]]][[.target.col[1]]] <- sapply(1:length(tc), function (k) {
        if (k > 1) {
          if (tc[k] == tc[k-1]) {
            paste0(tc[k], '.', dups.i[k])
          } else {
            tc[k]
          }
        } else {
          tc[k]
        }
      }, USE.NAMES = F)
    } else {
      dt.list[[names(.data)[i]]] <- do.call(data.table, sapply(c(.target.col, .col.name), function (cn) {
        x <- NA
        class(x) <- class(.data[[1]][1, cn])
        if (length(.target.col) == 1) {
          rep(x, times = length(.targets))
        } else {
          rep(x, times = nrow(.targets))
        }
      }, simplify = F))
      if (length(.target.col) == 1) {
        dt.list[[names(.data)[i]]][[.target.col]] <- .targets
      } else {
        for (tc in .target.col) {
          dt.list[[names(.data)[i]]][[tc]] <- .targets[[tc]]
        }
      }
      setnames(dt.list[[names(.data)[i]]], c(.target.col, .col.name))
      setkeyv(dt.list[[names(.data)[i]]], .target.col)
    }
  }
  if (.verbose) close(pb)
  res <- dt.list[[1]]
  if (length(dt.list) > 1) {
    for (i in 2:length(dt.list)) {
      res <- merge(res, dt.list[[i]], by = .target.col, all = T, allow.cartesian = T, suffixes = paste0('.', c(names(dt.list)[i-1], names(dt.list)[i])))
    }
  }
  for (i in 1:length(.col.name)) {
    if (.col.name[i] %in% names(res)) {
      setnames(res, .col.name[i], paste0(.col.name[i], '.', names(dt.list)[length(dt.list)]))
    }
  }
  res <- as.data.frame(res, stringsAsFactors = F)
  if (length(.col.name) > 1) {
    res <- res[, c(1:length(.target.col), cbind(seq(length(.target.col) + 1, ncol(res), 2), seq(length(.target.col) + 2, ncol(res), 2)))]
  }
  res[[.target.col[1]]] <- sapply(strsplit(res[[.target.col[1]]], '.', T, F, T), '[', 1, USE.NAMES = F)
  res <- res[order(res[[.target.col[1]]]),]
  
  tc <- res[[.target.col[1]]]
  dups <- duplicated(tc)
  dups.i <- rep('', length(dups))
  k <- 1
  for (j in 1:length(dups)) {
    if (dups[j]) {
      dups.i[j] <- paste0('.', as.character(k))
      k <- k + 1
    } else {
      k <- 1
    }
  }
  row.names(res) <- paste0(res[[.target.col[1]]], dups.i)
  
  res
}


#' Perform sequential cross starting from the top of a data frame.
#' 
#' @aliases top.cross top.cross.vec top.cross.plot
#' 
#' @description
#' \code{top.cross} - get top crosses of the given type between each pair of the given data.frames with \code{top.cross} function.
#' 
#' \code{top.cross.vec} - get vector of cross values for each top with the \code{top.cross.vec} function.
#' 
#' \code{top.cross.plot} - plot a plots with result with the \code{top.cross.plot} function.
#' 
#' @usage
#' top.cross(.data, .n = NA, .data2 = NULL, .type = 'ave', .norm = F, .verbose = T)
#' 
#' top.cross.vec(.top.cross.res, .i, .j)
#' 
#' top.cross.plot(.top.cross.res, .xlab = 'Top X clonotypes', 
#'                .ylab = 'Normalised number of shared clonotypes', .nrow = 2,
#'                .legend.ncol = 1, .logx = T, .logy = T)
#' 
#' @param .data Either list of data.frames or a data.frame.
#' @param .n Integer vector of parameter appled to the head function; same as .n in the top.fun function. See "Details" for more information.
#' @param .data2 Second data.frame or NULL if .data is a list.
#' @param .type Parameter .type to the \code{tcR::intersect} function.
#' @param .norm Parameter .norm to the \code{tcR::intersect} function.
#' @param .top.cross.res Result from the \code{top.cross} function.
#' @param .i,.j Coordinate of a cell in each matrix.
#' @param .xlab Name for a x-lab.
#' @param .ylab Name for a y-lab.
#' @param .nrow Number of rows of sub-plots in the output plot.
#' @param .legend.ncol Number of columns in the output legend.
#' @param .logx if T then transform x-axis to log-scale.
#' @param .logy if T then transform y-axis to log-scale.
#' @param .verbose if T then plot a progress bar.
#' 
#' @return
#' \code{top.cross} - return list for each element in \code{.n} with intersection matrix (from \code{tcR::intersectClonesets}).
#' 
#' \code{top.cross.vec} - vector of length \code{.n} with \code{.i}:\code{.j} elements of each matrix.
#' 
#' \code{top.cross.plot} - grid / ggplot object.
#' 
#' @details
#' Parameter \code{.n} can have two possible values. It could be either integer vector of numbers (same as in the \code{top.fun} function) or
#' NA and then it will be replaced internally by the value \code{.n <- seq(5000, min(sapply(.data, nrow)), 5000)}.
#' 
#' @seealso \code{\link{intersect}}
#' 
#' @examples
#' \dontrun{
#' immdata.top <- top.cross(immdata)
#' top.cross.plot(immdata.top)
#' }
top.cross <- function (.data, .n = NA, .data2 = NULL, .type = 'ave', .norm = F, .verbose = T) {
  if (!is.null(.data2)) { .data <- list(.data, .data2) }
  
  if (is.na(.n)[1]) { .n <- seq(5000, min(sapply(.data, nrow)), 5000) }
  
#   res <- lapply(.n, function(i) apply.symm(.data, cross, .head = i, .type = .type, .norm = .norm, .verbose=F))
  res <- lapply(.n, function(i) {
    if (.verbose) cat('Head == ', i, ' :\n', sep = '')
    intersectClonesets(.data, .type, .head = i, .norm = .norm, .verbose = .verbose)})
  names(res) <- .n
  res
}

top.cross.vec <- function (.top.cross.res, .i, .j) {
  sapply(.top.cross.res, function (mat) mat[.i, .j] )
}

top.cross.plot <- function (.top.cross.res, .xlab = 'Top X clonotypes', .ylab = 'Normalised number of shared clonotypes', .nrow = 2, .legend.ncol = 1, .logx = T, .logy = T) {
  data.names <- colnames(.top.cross.res[[1]])
  len <- length(data.names)
  ps <- lapply(1:len, function (i) { 
    data <- data.frame(Top = as.numeric(names(.top.cross.res)), sapply((1:len), function (j) {
      res <- top.cross.vec(.top.cross.res, i, j)
#       res[is.na(res)] <- 0
      res
    }), sapply((1:len), function (j) {
      rep(data.names[j], times=length(.top.cross.res))
    }))
    names(data) <- c('Top', data.names, paste0(data.names, 'C'))
    p <- ggplot(data = data)
    for (j in (1:len)) {
      p <- p + 
        geom_point(aes_string(x = 'Top', y = data.names[j], color = paste0(data.names[j], 'C'))) +
        geom_line(aes_string(x = 'Top', y = data.names[j], color = paste0(data.names[j], 'C'))) +
        xlab(.xlab) + ylab(.ylab) + theme_linedraw()
    }
    p + ggtitle(data.names[i]) + theme(legend.position = 'none')
    if (.logx) {
      p <- p + scale_x_log10()
    }
    if (.logy) {
      p <- p + scale_y_log10()
    }
    p <- p + .colourblind.discrete(len, T)
  } )  
  
  data <- data.frame(Top = as.numeric(names(.top.cross.res)), 
                     sapply((1:len), function (j) rep.int(1, length(.top.cross.res))), 
                     sapply((1:len), function (j) rep(data.names[j], times=length(.top.cross.res))))
  names(data) <- c('Top', data.names, paste0(data.names, 'C'))
  p <- ggplot(data = data)
  for (j in (1:len)) {
    p <- p + geom_point(aes_string(x = 'Top', y = data.names[j], color = paste0(data.names[j], 'C')))
  }
  
  #   add.legend(ps, .nrow, .legend.ncol)
  sample.plot <- p + .colourblind.discrete(len, T)
  
  leg <- gtable_filter(ggplot_gtable(ggplot_build(sample.plot + guides(colour=guide_legend(title = 'Subject', ncol=.legend.ncol)))), "guide-box")
  grid.arrange(do.call(arrangeGrob, c(lapply(ps, function (x) x + theme(legend.position="none")), nrow = .nrow)), leg, widths=unit.c(unit(1, "npc") - leg$width, leg$width), nrow = 1, top ='Top crosses')
#   arrangeGrob(do.call(arrangeGrob, c(ps[1:2], nrow = .nrow)), leg, widths=unit.c(unit(1, "npc") - leg$width, leg$width), nrow = 1, top ='Top cross')
}


#' Bootstrap for data frames in package tcR.
#' 
#' @description
#' Resample rows (i.e., clones) in the given data frame and apply the given function to them.
#' 
#' @param .data Data frame.
#' @param .fun Function applied to each sample.
#' @param .n Number of iterations (i.e., size of a resulting distribution).
#' @param .size Size of samples. For \code{.sim} == "uniform" stands for number of rows to take.
#' For \code{.sim} == "percentage" stands for number of UMIs / read counts to take.
#' @param .sim A character string indicating the type of simulation required. Possible values are
#' "uniform" or "percentage". See "Details" for more details of type of simulation.
#' @param .postfun Function applied to the resulting list: list of results from each processed sample.
#' @param .verbose if T then show progress bar.
#' @param .prop.col Column with proportions for each clonotype.
#' @param ... Further values passed to \code{.fun}.
#' 
#' @return
#' Either result from \code{.postfun} or list of length \code{.n} with values of \code{.fun}.
#' 
#' @details
#' Argument \code{.sim} can take two possible values: "uniform" (for uniform distribution), when
#' each row can be taken with equal probability, and "perccentage" when each row can be taken with
#' probability equal to its "Read.proportion" column.
#' 
#' @examples
#' \dontrun{
#' # Apply entropy.seg function to samples of size 20000 from immdata$B data frame for 100 iterations.
#' bootstrap.tcr(immdata[[2]], .fun = entropy.seg, .n = 100, .size = 20000, .sim = 'uniform')
#' }
bootstrap.tcr <- function (.data, .fun = entropy.seg, .n = 1000,
                           .size = nrow(.data), .sim = c('uniform', 'percentage'),
                           .postfun = function (x) { unlist(x) }, .verbose = T, .prop.col = 'Read.proportion',
                           ...) {
  
  .sample.fun <- function (d) {
    d[sample(1:nrow(d), .size, T),]
  }
  
  if (.sim[1] == 'percentage') {
    .sample.fun <- function (d) {
      new.reads <- rmultinom(1, .size, d[, .prop.col])
      d$Read.count <- new.reads
      d[, .prop.col] <- new.reads / sum(new.reads)
      d[new.reads > 0,]
    }
  }
  
  if (.verbose) { pb <- set.pb(.n) }
  res <- lapply(1:.n, function (n) {
    if (.verbose) { add.pb(pb) }
    .fun(.sample.fun(.data), ...)
  })
  if (.verbose) { close(pb) }

  .postfun(res)
}


# mean + IQR
#' Clonal space homeostasis.
#' 
#' @description
#' Compute clonal space homeostatsis - statistics of how many space occupied by clones
#' with specific proportions.
#' 
#' @param .data Cloneset data frame or list with such data frames.
#' @param .clone.types Named numeric vector.
#' @param .prop.col Which column to use for counting proportions.
#' 
#' @seealso \link{vis.clonal.space}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' # Compute summary space of clones, that occupy
#' # [0, .05) and [.05, 1] proportion.
#' clonal.space.homeostasis(twb, c(Low = .05, High = 1)))
#' #        Low (0 < X <= 0.05) High (0.05 < X <= 1)
#' # Subj.A           0.9421980           0.05780198
#' # Subj.B           0.9239454           0.07605463
#' # Subj.C           0.8279270           0.17207296
#' # Subj.D           1.0000000           0.00000000
#' # I.e., for Subj.D sum of all read proportions for clones
#' # which have read proportion between 0 and .05 is equal to 1.
#' }
clonal.space.homeostasis <- function (.data, .clone.types = c(Rare = .00001,
                                                              Small = .0001,
                                                              Medium = .001,
                                                              Large = .01,
                                                              Hyperexpanded = 1),
                                      .prop.col = 'Read.proportion') {
  .clone.types <- c(None = 0, .clone.types)
  
  if (has.class(.data, 'data.frame')) {
    .data <- list(Sample = .data)
  }
  
  mat <- matrix(0, length(.data), length(.clone.types) - 1, dimnames = list(names(.data), names(.clone.types)[-1]))
  .data <- lapply(.data, '[[', .prop.col)
  for (i in 2:length(.clone.types)) {
    mat[,i-1] <- sapply(.data, function (x) sum(x[x > .clone.types[i-1] & x <= .clone.types[i]]))
    colnames(mat)[i-1] <- paste0(names(.clone.types[i]), ' (', .clone.types[i-1], ' < X <= ', .clone.types[i], ')')
  }
  
  mat
}