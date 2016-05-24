########## Intersections among sets of sequences ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Resamp.values", "Fun.value", "P.value"))
}


#' Intersection between sets of sequences or any elements.
#' 
#' @aliases intersectClonesets intersectCount intersectLogic intersectIndices
#' 
#' @description
#' Functions for the intersection of data frames with TCR / Ig data. 
#' See the \code{repOverlap} function for a general interface to all overlap analysis functions.
#' 
#' \code{intersectClonesets} - returns number of similar elements in the given two clonesets / data frames or matrix
#' with counts of similar elements among each pair of objects in the given list.
#' 
#' \code{intersectCount} - similar to \code{tcR::intersectClonesets}, but with fewer parameters and only for two objects.
#' 
#' \code{intersectIndices} - returns matrix M with two columns, where element with index M[i, 1] in the first
#' given object is similar to an element with index M[i, 2] in the second given object.
#' 
#' \code{intersectLogic} - returns logic vector with TRUE values in positions, where element in the first given data frame
#' is found in the second given data frame.
#' 
#' @usage
#' intersectClonesets(.alpha = NULL, .beta = NULL, .type = "n0e", .head = -1, .norm = F,
#'           .verbose = F)
#' 
#' intersectCount(.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL)
#' 
#' intersectIndices(.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL)
#' 
#' intersectLogic(.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL)
#' 
#' @param .alpha Either first vector or data.frame or list with data.frames.
#' @param .beta Second vector or data.frame or type of intersection procedure (see the \code{.type} parameter) if \code{.alpha} is a list.
#' @param .type Types of intersection procedure if \code{.alpha} and \code{.beta} is data frames. String with 3 characters (see 'Details' for more information).
#' @param .head Parameter for the \code{head} function, applied before intersecting.
#' @param .method Method to use for intersecting string elements: 'exact' for exact matching, 'hamm' for matching strings which have <= 1 hamming distance,
#' 'lev' for matching strings which have <= 1 levenshtein (edit) distance between them.
#' @param .col Which columns use for fetching values to intersect. First supplied column matched with \code{.method}, others as exact values.
#' @param .norm If TRUE than normalise result by product of length or nrows of the given data.
#' @param .verbose if T then produce output of processing the data.
#' 
#' @details
#' Parameter \code{.type} of the \code{intersectClonesets} function is a string of length 3
#' [0an][0vja][ehl], where:
#' \enumerate{
#'  \item First character defines which elements intersect ("a" for elements from the column "CDR3.amino.acid.sequence", 
#'  "n" for elements from the column "CDR3.nucleotide.sequence", other characters - intersect elements as specified);
#'  \item Second character defines which columns additionaly script should use
#' ('0' for cross with no additional columns, 'v' for cross using the "V.gene" column, 
#' 'j' for cross using "J.gene" column, 'a' for cross using both "V.gene" and "J.gene" columns);
#'  \item Third character defines a method of search for similar sequences is use:
#'  "e" stands for the exact match of sequnces, "h" for match elements which have the Hamming distance between them
#'  equal to or less than 1, "l" for match elements which have the Levenshtein distance between tham equal to or less than 1.
#' }
#' 
#' @seealso  \link{repOverlap}, \link{vis.heatmap}, \link{ozScore}, \link{permutDistTest}, \link{vis.group.boxplot}
#' 
#' @return
#' \code{intersectClonesets} returns (normalised) number of similar elements or matrix with numbers of elements.
#' 
#' \code{intersectCount} returns number of similar elements.
#' 
#' \code{intersectIndices} returns 2-row matrix with the first column stands for an index of an element in the given \code{x}, and the second column stands for an index of an element of \code{y} which is similar to a relative element in \code{x}; 
#' 
#' \code{intersectLogic} returns logical vector of \code{length(x)} or \code{nrow(x)}, where TRUE at position \code{i} means that element with index {i} has been found in the \code{y}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' # Equivalent to intersectClonesets(twb[[1]]$CDR3.nucleotide.sequence,
#' #                         twb[[2]]$CDR3.nucleotide.sequence)
#' # or intersectCount(twb[[1]]$CDR3.nucleotide.sequence,
#' #                    twb[[2]]$CDR3.nucleotide.sequence)
#' # First "n" stands for a "CDR3.nucleotide.sequence" column, "e" for exact match.
#' twb.12.n0e <- intersectClonesets(twb[[1]], twb[[2]], 'n0e')
#' stopifnot(twb.12.n0e == 46)
#' # First "a" stands for "CDR3.amino.acid.sequence" column.
#' # Second "v" means that intersect should also use the "V.gene" column.
#' intersectClonesets(twb[[1]], twb[[2]], 'ave')
#' # Works also on lists, performs all possible pairwise intersections.
#' intersectClonesets(twb, 'ave')
#' # Plot results.
#' vis.heatmap(intersectClonesets(twb, 'ave'), .title = 'twb - (ave)-intersection', .labs = '')
#' # Get elements which are in both twb[[1]] and twb[[2]].
#' # Elements are tuples of CDR3 nucleotide sequence and corresponding V-segment
#' imm.1.2 <- intersectLogic(twb[[1]], twb[[2]],
#'                            .col = c('CDR3.amino.acid.sequence', 'V.gene'))  
#' head(twb[[1]][imm.1.2, c('CDR3.amino.acid.sequence', 'V.gene')])
#' data(twb)
#' ov <- repOverlap(twb)
#' sb <- matrixSubgroups(ov, list(tw1 = c('Subj.A', 'Subj.B'), tw2 = c('Subj.C', 'Subj.D')));
#' vis.group.boxplot(sb)
#' }
intersectClonesets <- function (.alpha = NULL, .beta = NULL, .type = 'n0e', .head = -1, .norm = F, .verbose = F) {
  if (class(.alpha) == 'list') {
    if (class(.beta) == 'character') {
      .type <- .beta
    }
    apply.symm(.alpha, intersectClonesets, .head = .head, .type = .type, .norm = .norm, .verbose = .verbose)
  } else {
    if (.head != -1) {
      .alpha <- head(.alpha, .head)
      .beta <- head(.beta, .head)
    }
    
    cols <- NULL
    if (class(.alpha) == 'data.frame') {
      if (substr(.type, 1, 1) == 'a' || substr(.type, 1, 1) == 'n') {
        if (substr(.type, 1, 1) == 'a') {
          cols <- 'CDR3.amino.acid.sequence'
        } else {
          cols <- 'CDR3.nucleotide.sequence'
        }
        
        if (substr(.type, 2, 2) == 'v') {
          cols <- c(cols, 'V.gene')
        } else if (substr(.type, 2, 2) == 'j') {
          cols <- c(cols, 'J.gene')
        } else if (substr(.type, 2, 2) == 'a') {
          cols <- c(cols, 'V.gene', 'J.gene')
        } else if (substr(.type, 2, 2) != '0') {
          cat("Second character in .type:", .type, 'is unknown!\n')
        }
      }
    }
    
    method <- switch(substr(.type, 3, 3),
                  e = 'exact',
                  h = 'hamm',
                  l = 'lev')
    res <- intersectCount(.alpha, .beta, method, cols)
    if (.norm) {
      if (is.null(dim(.alpha))) {
        res <- res / (as.numeric(length(.alpha)) * length(.beta))
        # res <- res * max(length(.alpha), length(.beta))
      } else {
        res <- res / (as.numeric(nrow(.alpha)) * nrow(.beta))
        # res <- res * max(nrow(.alpha), nrow(.beta))
      }
    }
    res
  }
}

intersectCount <- function (.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL) {
  if (.method[1] == 'exact' && (is.null(.col) || length(.col) == 1)) {
    if (is.null(.col)) {
      length(base::intersect(.alpha, .beta))
    } else {
      length(base::intersect(.alpha[,.col], .beta[,.col]))
    }
  } else {
    if (.method[1] == 'exact') {
      nrow(dplyr::intersect(.alpha[, .col], .beta[, .col]))
    } else {
      res <- intersectIndices(unique(.alpha), unique(.beta), .method, .col)
      if (!is.na(res[1,1])) {
        nrow(res)
      } else {
        0
      }
    }
  }
}

intersectIndices <- function (.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL) {    
  if (!is.null(.col)) {
    .alpha <- as.data.frame(.alpha, stringsAsFactors = F)
    .beta <- as.data.frame(.beta, stringsAsFactors = F)
    if (length(.col) == 1) {
      .alpha <- .alpha[, .col]
      .beta <- .beta[, .col]
    } else {
      if (length(.col) == 2) {
        alpha.cols <- as.matrix(.alpha[, .col[2]])
        beta.cols <- as.matrix(.beta[, .col[2]])
      } else {
        alpha.cols <- .alpha[, .col[2:length(.col)]]
        beta.cols <- .beta[, .col[2:length(.col)]]
      }
      .alpha <- .alpha[, .col[1]]
      .beta <- .beta[, .col[1]]
    }
  }
  
  res <- find.similar.sequences(.alpha, .beta, .method, 1, F)
  
  if (is.na(res[1,1])) {
    return(res)
  }
  
  if (!is.null(.col) && length(.col) > 1) {
    for (i in 1:ncol(alpha.cols)) {
      res <- res[alpha.cols[res[,1], i] == beta.cols[res[,2], i], ]
    }
  }
  
  if(!is.matrix(res)) {
    res <- rbind(res)
  }
  
  if (dim(res)[1] == 0) {
    res <- matrix(c(NA, NA), 1, 2)
  }
  
  res
}

intersectLogic <- function (.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL) {
  ind <- intersectIndices(.alpha, .beta, .method, .col)
  if (is.null(.col)) {
    logic <- rep(FALSE, length(.alpha))
  } else {
    logic <- rep(FALSE, nrow(.alpha))
  }
  logic[ind[,1]] <- TRUE
  logic
}


#' Compute convergence characteristics of repertoires.
#' 
#' @description 
#' Get a number of rows with similar aminoacid sequence but different nucleotide sequence.
#' 
#' @param .alpha Either data frame with columns \code{.col.nuc} and \code{.col.aa} or list with such data frames.
#' @param .beta Either data frame or none.
#' @param .col.nuc Name of the column with nucleotide sequences.
#' @param .col.aa Name of the columnw ith aminoacid sequences.
#' @return If \code{.alpha} is data frame, than integer vector of length 2 with . If \code{.alpha} is a list
#' than matrix M with M[i,j] = convergence.index(.alpha[[i]], .alpha[[j]]).
convergence.index <- function (.alpha, .beta, .col.nuc = 'CDR3.nucleotide.sequence', .col.aa = 'CDR3.amino.acid.sequence') {
  if (has.class(.alpha, 'list')) {
    return(apply.asymm(.alpha, function (x,y) convergence.index(x, y)[1]))
  }
  
  a.nuc.logic <- .alpha[, .col.nuc] %in% .beta[, .col.nuc]
  b.nuc.logic <- .beta[, .col.nuc] %in% .alpha[, .col.nuc]
  a.aa.logic <- .alpha[, .col.aa] %in% .beta[, .col.aa]
  b.aa.logic <- .beta[, .col.aa] %in% .alpha[, .col.aa]
  c(Amino.acid.not.nucleotide.count.1 = length(unique(.alpha[a.aa.logic & !(a.nuc.logic), .col.aa])),
    Amino.acid.not.nucleotide.count.2 = length(unique(.beta[b.aa.logic & !(b.nuc.logic), .col.aa])))
}


#' Overlap Z-score.
#' 
#' @description
#' Compute OZ-scores ("overlap Z scores") for values in the given matrix of overlaps, i.e.,.
#' for each value compute the number of standart deviations from the mean of the matrix.
#' 
#' @param .mat Matrix with overlap values.
#' @param .symm If T then remove lower triangle matrix from counting. Doesn't work if the matrix
#' has different number of rows and columns.
#' @param .as.matrix If T then return
#' @param .val.col If .as.matrix is T then this is a name of the column to build matrix upon:
#' either "oz" for the OZ-score column, "abs" for the absolute OZ-score column, or "norm" for the
#' normalised absolute OZ-score column.
#' 
#' @seealso \link{repOverlap}, \link{intersectClonesets}, \link{permutDistTest}
#' 
#' @examples 
#' \dontrun{
#' data(twb)
#' mat <- repOverlap(twb)
#' ozScore(mat)
#' # Take 3x3 matrix
#' ozScore(mat[1:3, 1:3])
#' # Return as matrix with OZ scores
#' ozmat <- ozScore(mat, T, T, 'oz')
#' # Return as matrix with normalised absolute OZ scores
#' oznorm <- ozScore(mat, T, T, 'norm')
#' # Plot it as boxplots
#' sb <- matrixSubgroups(oznorm, list(tw1 = c('Subj.A', 'Subj.B'), tw2 = c('Subj.C', 'Subj.D')));
#' vis.group.boxplot(sb)
#' }
ozScore <- function (.mat, .symm = T, .as.matrix = F, .val.col = c('norm', 'abs', 'oz')) {
  if (.symm && nrow(.mat) == ncol(.mat)) {
    .mat[lower.tri(.mat)] <- NA
  }
  
  res <- melt(.mat, na.rm = T)
  res$OZ <- NA
  colnames(res) <- c("Rep1", "Rep2", "Overlap", "OZ")
  
  res$OZ <- (res$Overlap - mean(res$Overlap)) / sd(res$Overlap)
  res$Abs.OZ <- abs(res$OZ)
  res$Norm.abs.OZ <- (res$Abs.OZ - min(res$Abs.OZ)) / (max(res$Abs.OZ) - min(res$Abs.OZ))
  res <- res[order(res$Abs.OZ, decreasing = T),]
  row.names(res) <- NULL
  
  if (.as.matrix) {
    res <- acast(res, Rep1 ~ Rep2, value.var = switch(.val.col[1], oz = 'OZ', abs = 'Abs.OZ', norm = 'Norm.abs.OZ'))
    if (sum(!colnames(.mat) %in% colnames(res))) {
      res <- cbind(NA, res)
      colnames(res)[1] <- colnames(.mat)[!(colnames(.mat) %in% colnames(res))]
    }
    
    if (sum(!row.names(.mat) %in% row.names(res))) {
      res <- rbind(res, NA)
      row.names(res)[nrow(res)] <- row.names(.mat)[!(row.names(.mat) %in% row.names(res))]
    }
  }
  
  res
}


#' Monte Carlo permutation test for pairwise and one-vs-all-wise within- and inter-group differences in a set of repertoires.
#' 
#' @description
#' WARNING: this is an experimental procedure, work is still in progress.
#' 
#' Perform permutation tests of distances among groups for the given groups of samples and matrix of distances among all samples.
#' 
#' @param .mat Symmetric matrix of repertoire distances.
#' @param .groups Named list with names of repertoires in groups.
#' @param .n Number of permutations for each pair of group.
#' @param .fun A function to apply to distances.
#' @param .signif Significance level. Below this value hypotheses counts as significant.
#' @param .plot If T than plot the output results. Else return them as a data frame.
#' @param .xlab X lab label.
#' @param .title Main title of the plot.
#' @param .hjust Value for adjusting the x coordinate of p-value labels on plots.
#' @param .vjust Value for adjusting the y coordinate of p-value labels on plots.
#' 
#' @seealso \link{repOverlap}, \link{intersectClonesets}, \link{ozScore}, \link{pca2euclid}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' mat <- repOverlap(twb)
#' permutDistTest(mat, list(tw1 = c('Subj.A', 'Subj.B'), tw2 = c('Subj.C', 'Subj.D')))
#' permutDistTest(mat, list(tw1 = c('Subj.A', 'Subj.B'), tw2 = c('Subj.C', 'Subj.D')), .fun = median)
#' }
permutDistTest <- function (.mat, .groups, .n = 1000, .fun = mean, .signif = .05,
                           .plot = T, .xlab = "Values", .title = "Monte Carlo permutation testing of overlaps",
                           .hjust = -.1, .vjust = -4) {
  
  cat("WARNING: this is an experimental procedure, work is still in progress.")
  
  .pairwise.test <- function (.mat, .group.logic, .n, .fun) {
    within.val.gr1 = .fun(.mat[.group.logic, .group.logic][upper.tri(.mat[.group.logic, .group.logic])])
    within.val.gr2 = .fun(.mat[!.group.logic, !.group.logic][upper.tri(.mat[!.group.logic, !.group.logic])])
    inter.val = .fun(.mat[.group.logic, !.group.logic])
    
    within.resamp.gr1 = c()
    within.resamp.gr2 = c()
    inter.resamp = c()
    
    gr1.size <- sum(.group.logic)
    gr2.size <- sum(!.group.logic)
    
    for (iter in 1:.n) {
      new.group.logic <- rep(F, nrow(.mat))
      new.group.logic[ sample(1:nrow(.mat), gr1.size, F) ] <- T
      
      within.resamp.gr1 = c(within.resamp.gr1, .fun(.mat[new.group.logic, new.group.logic][upper.tri(.mat[new.group.logic, new.group.logic])]))
      within.resamp.gr2 = c(within.resamp.gr2, .fun(.mat[!new.group.logic, !new.group.logic][upper.tri(.mat[!new.group.logic, !new.group.logic])]))
      inter.resamp = c(inter.resamp, .fun(.mat[new.group.logic, !new.group.logic]))
    }
    
    list(within.val.gr1 = within.val.gr1,
         within.val.gr2 = within.val.gr2,
         inter.val = inter.val,
         within.resamp.gr1 = within.resamp.gr1,
         within.resamp.gr2 = within.resamp.gr2,
         inter.resamp = inter.resamp,
         within.p.gr1 = sum(within.val.gr1 > within.resamp.gr1) / .n,
         within.p.gr2 = sum(within.val.gr2 > within.resamp.gr2) / .n,
         inter.p = sum(inter.val > inter.resamp) / .n)
  }
  
  .intergroup.name <- function (.gr1, .gr2) {
    paste0(.gr1, ':', .gr2)
  }
  
  res <- list()
  
  pb <- set.pb((length(.groups) - 1) * length(.groups) / 2)
  
  k <- 1
  
  p.vals <- list()
  for (gr in names(.groups)) {
    p.vals[[gr]] <- list(within = c(), inter = c())
  }
  
  for (i in 1:(length(.groups)-1)) {
    gr1 <- names(.groups)[i]
    for (j in (i+1):length(.groups)) {
      gr2 <- names(.groups)[j]
      
      submat <- .mat[ row.names(.mat) %in% .groups[[i]] | row.names(.mat) %in% .groups[[j]] , colnames(.mat) %in% .groups[[i]] | colnames(.mat) %in% .groups[[j]] ]
      submat <- submat[, row.names(submat)]
      group.logic <- row.names(submat) %in% .groups[[i]]
      
      test.res.list <- .pairwise.test(submat, group.logic, .n, .fun)
      
      res[[k]] <- rbind(data.frame(Group1 = gr1,
                                   Group2 = gr1,
                                   Fun.value = test.res.list$within.val.gr1, 
                                   Resamp.values = test.res.list$within.resamp.gr1,
                                   P.value = test.res.list$within.p.gr1,
                                   stringsAsFactors = F), 
                        data.frame(Group1 = gr1,
                                   Group2 = gr2,
                                   Fun.value = test.res.list$inter.val, 
                                   Resamp.values = test.res.list$inter.resamp, 
                                   P.value = test.res.list$inter.p,
                                   stringsAsFactors = F), 
                        data.frame(Group1 = gr2,
                                   Group2 = gr2,
                                   Fun.value = test.res.list$within.val.gr2, 
                                   Resamp.values = test.res.list$within.resamp.gr2, 
                                   P.value = test.res.list$within.p.gr2,
                                   stringsAsFactors = F))
      k <- k + 1
      
      p.vals[[gr1]][["within"]] <- c(p.vals[[gr1]][["within"]], test.res.list$within.p.gr1)
      names(p.vals[[gr1]][["within"]])[length(p.vals[[gr1]][["within"]])] <- gr2
      
      p.vals[[gr2]][["within"]] <- c(p.vals[[gr2]][["within"]], test.res.list$within.p.gr2)
      names(p.vals[[gr2]][["within"]])[length(p.vals[[gr2]][["within"]])] <- gr1
      
      p.vals[[gr1]][["inter"]] <- c(p.vals[[gr1]][["inter"]], test.res.list$inter.p)
      names(p.vals[[gr1]][["inter"]])[length(p.vals[[gr1]][["inter"]])] <- gr2
      
      p.vals[[gr2]][["inter"]] <- c(p.vals[[gr2]][["inter"]], test.res.list$inter.p)
      names(p.vals[[gr2]][["inter"]])[length(p.vals[[gr2]][["inter"]])] <- gr1
      
      add.pb(pb)
    }
  }
  close(pb)
  
  p <- list()
  for (i in 1:(k-1)) {
    faclevels <- unique(c(res[[i]]$Group1, res[[i]]$Group2))
    faclevels <- c(faclevels[1], paste0(faclevels[1], " ~ ", faclevels[2]), faclevels[2])
    
    res[[i]]$Group3 <- paste0(res[[i]]$Group1, " ~ ", res[[i]]$Group2)
    res[[i]]$Group3[res[[i]]$Group1 == res[[i]]$Group2] <- res[[i]]$Group1[res[[i]]$Group1 == res[[i]]$Group2]
    res[[i]]$Group3 <- factor(res[[i]]$Group3, faclevels)
    
    tmp <- do.call(rbind, lapply(split(res[[i]], res[[i]]$Group3), function (x) x[1,]))
    
    p[[i]] <- ggplot(res[[i]], aes(x = Resamp.values)) + 
      geom_vline(aes(xintercept = Fun.value), size = 1, linetype = "twodash", colour = "blue") + 
      geom_histogram(alpha = .4, colour = 'grey40') + 
      geom_text(aes(x = 0, y = 0, label = paste0("P = ", P.value)), hjust = .hjust, vjust = .vjust, size = 5, data = tmp) +
      facet_wrap(~ Group3, ncol = 3) + 
      ylab("Permutations, num") + 
      xlab(.xlab) +
      theme_bw()
  }
  
  fmt.width <- 0
  for (gr.name in names(.groups)) {
    fmt.width <- max(fmt.width, nchar(gr.name))
  }
  fmt.width <- fmt.width + 2
  
  flag <- F
  cat("Significant differences found (OBS - observed, SIM - simulated):\n")
  for (gr.name in names(p.vals)) {
    
    w.signs <- c()
    i.signs <- c()
    
    for (i in 1:length(p.vals[[gr.name]][["within"]])) {
      if (p.vals[[gr.name]][["within"]][i] < 1 - p.vals[[gr.name]][["within"]][i]) {
        w.signs <- c(w.signs, ">")
      } else {
        w.signs <- c(w.signs, "<")
        p.vals[[gr.name]][["within"]][i] <- 1 - p.vals[[gr.name]][["within"]][i]
      }
    }
    for (i in 1:length(p.vals[[gr.name]][["inter"]])) {
      if (p.vals[[gr.name]][["inter"]][i] < 1 - p.vals[[gr.name]][["inter"]][i]) {
        i.signs <- c(i.signs, ">")
      } else {
        i.signs <- c(i.signs, "<")
        p.vals[[gr.name]][["inter"]][i] <- 1 - p.vals[[gr.name]][["inter"]][i]
      }
    }
    
    p.vals[[gr.name]][["within"]] <- p.adjust(p.vals[[gr.name]][["within"]], "BH")
    p.vals[[gr.name]][["inter"]] <- p.adjust(p.vals[[gr.name]][["inter"]], "BH")
    
#     cat("\tThis groups may have come from different populations\n")
#     cat("\tThis groups may have non-random overlap and come from different populations\n")
    
    for (i in 1:length(p.vals[[gr.name]][["within"]])) {
      pval <- p.vals[[gr.name]][["within"]][i]
      if (pval <= .signif) {
        cat('  Within  ', 
            formatC(paste0('"', gr.name, '"'), width = -fmt.width), 
            ' in a pool with ', 
            formatC(paste0('"', names(p.vals[[gr.name]][["within"]])[i], '"'), width = fmt.width), 
            '  :  P(OBS ', w.signs[i], ' SIM) = ', pval, "\n", sep ='')
        flag <- T
      }
    }
    
    for (i in 1:length(p.vals[[gr.name]][["inter"]])) {
      pval <- p.vals[[gr.name]][["inter"]][i]
      if (pval <= .signif) {
        if (pval <= .signif) {
          cat('  Between ', 
              formatC(paste0('"', gr.name, '"'), width = -fmt.width), 
              ' and ', 
              formatC(paste0('"', names(p.vals[[gr.name]][["inter"]])[i], '"'), width = fmt.width + 11), 
              '  :  P(OBS ', i.signs[i], ' SIM) = ', pval, "\n", sep = '')
          flag <- T
        }
      }
    }
  }
  if (!flag) { cat("No significant differences have been found.\n") }
  
  if (.plot) {
    suppressMessages(do.call(grid.arrange, c(p, ncol = 1, top = .title)))
  } else {
    res
  }
}