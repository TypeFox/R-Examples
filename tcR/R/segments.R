########## Statistics and analysis of Variable and Joining genes usage ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c('PC1', 'PC2', "Subject"))
}


#' Gene usage.
#' 
#' @aliases geneUsage
#' 
#' @description 
#' Compute frequencies or counts of gene segments ("V / J - usage").
#' 
#' @param .data Cloneset data frame or a list with clonesets.
#' @param .genes Either one of the gene alphabet (e.g., HUMAN_TRBV, \link{genealphabets}) or list with two gene alphabets for computing 
#' joint distribution.
#' @param .quant Which column to use for the quantity of clonotypes: NA for computing only number of genes without using clonotype counts, 
#' "read.count" for the "Read.count" column, "umi.count" for the "Umi.count" column, "read.prop" for the "Read.proportion" column,
#' "umi.prop" for the "Umi.proportion" column.
#' @param .norm If T then return proportions of resulting counting of genes.
#' @param .ambig If F than remove from counting genes which are not presented in the given gene alphabet(s).
#' 
#' @return 
#' If \code{.data} is a cloneset and \code{.genes} is NOT a list than return a data frame with first column "Gene" with genes and second with counts / proportions.
#' 
#' If \code{.data} is a list with clonesets and \code{.genes} is NOT a list than return a data frame with first column "Gene" 
#' with genes and other columns with counts / proportions for each cloneset in the input list.
#' 
#' If \code{.data} is a cloneset and \code{.genes} IS a list than return a matrix with gene segments for the first gene in \code{.genes}
#' and column names for the second gene in \code{.genes}. See "Examples".
#' 
#' If \code{.data} is a list with clonesets and \code{.genes} IS a list than return a list with matrices like in the previous case.
#' 
#' @seealso \code{\link{genealphabets}}, \code{\link{vis.gene.usage}}, \code{\link{pca.segments}}
#' 
#' @examples
#' \dontrun{
#' # Load your data
#' data(twb)
#' # compute V-segments frequencies of human TCR beta.
#' seg <- geneUsage(twb, HUMAN_TRBV, .norm = T)
#' # plot V-segments frequencies as a heatmap
#' vis.heatmap(seg, .labs = c("Sample", "V gene"))
#' # plot V-segments frequencies directly from clonesets
#' vis.gene.usage(twb, HUMAN_TRBV)
#' # plot V-segments frequencies from the gene frequencies
#' vis.gene.usage(seg, NA)
#' # Compute V-J joint usage.
#' geneUsage(twb, list(HUMAN_TRBV, HUMAN_TRBJ))
#' # for future:
#' # geneUsage(twb, "human", "trbv")
#' }
geneUsage <- function (.data, .genes = HUMAN_TRBV_MITCR, .quant = c(NA, "read.count", "umi.count", "read.prop", "umi.prop"), 
                       .norm = F, .ambig = F #, .species = c("human", "mouse"), .genes = c("trbv", "trbd", "trbj")
                       ) {
  
  .process.df <- function (.df, .quant, .cols) {
    cast.fun <- dcast
    if (length(.cols) == 2) { cast.fun <- acast; len <- 2 }
    
    for (i in 1:length(.cols)) {
      .df[ which(!(.df[[.cols[i]]] %in% .genes[[i]])), .cols[i] ] <- "Ambiguous"
    }
    
    count.fun <- "n()"
    if (!is.na(.quant)) { count.fun <- paste0("sum(", .quant, ")", collapse = "", sep = "")}
    
    if (length(.cols) == 1) { .cols <- c(.cols, '.'); len <- 1}
    
    cast.fun(summarise_(grouped_df(select_(.df, .dots = as.list(na.exclude(c(.quant, .cols[1:len])))), lapply(.cols[1:len], as.name)), Freq = count.fun), as.formula(paste0(.cols[1], " ~ ", .cols[2])), value.var = 'Freq')
  }
  
  
  .fix.ambig <- function (.res, .ambig) {
    if (length(.genes) == 2) {
      .res <- lapply(.res, function (x) {
        x[row.names(x) != "Ambiguous", ][, colnames(x) != "Ambiguous"]
      })
      if (length(.data) == 1) {
        .res <- .res[[1]]
      }
      .res
    } else {
      .res[.res[,1] != "Ambiguous", ]
    }
  }
  
  
  quant <- NA
  if (!is.na(.quant[1])) { quant <- .column.choice(.quant, T) }
  
  if (has.class(.data, 'data.frame')) { .data <- list(Sample = .data) }
  
  if (has.class(.genes, 'list')) {    
    genecols <- c(paste0(substr(.genes[[1]][1], 4, 4), ".gene"), paste0(substr(.genes[[2]][1], 4, 4), ".gene"))
  } else {
    genecols <- paste0(substr(.genes[1], 4, 4), ".gene")
    .genes <- list(.genes)
  }
  
  tbls <- lapply(.data, .process.df, .quant = quant, .cols = genecols)
  
  
  # JOINT GENE DISTRIBUTION
  if (length(.genes) == 2) {
    tbls <- lapply(tbls, function (x) {
      genrows <- .genes[[1]][is.na(match(.genes[[1]], row.names(x)))]
      gencols <- .genes[[2]][is.na(match(.genes[[2]], colnames(x)))]
      if (length(genrows) > 0) {
        x <- do.call(rbind, c(list(x), rep.int(0, length(genrows))))
        row.names(x)[(nrow(x) - length(genrows) + 1):nrow(x)] <- genrows
      }
      
      if (length(gencols) > 0) {
        x <- do.call(cbind, c(list(x), rep.int(0, length(gencols))))
        colnames(x)[(ncol(x) - length(gencols) + 1):ncol(x)] <- gencols
      }

      x[is.na(x)] <- 0
      x[order(.genes[[1]]), ][, order(.genes[[2]])]
    })
    
    if (.norm) {
      tbls <- lapply(tbls, function (x) x / sum(x))
    }
    return(.fix.ambig(tbls, .ambig))
  }
  
  # SINGLE GENE DISTRIBUTION
  res <- tbls[[1]]
  colnames(res) <- c("Gene", names(.data)[1])
  
  if (length(.data) > 1) {
    for (i in 2:length(.data)) {
      colnames(tbls[[i]]) <- c("Gene", names(.data)[i])
      res <- merge(res, tbls[[i]], by = "Gene", all = T)
    }
  }
  res <- merge(res, data.frame(Gene = .genes[[1]], Something = 0, stringsAsFactors = F), by = "Gene", all = T)
  res <- res[, -ncol(res)]
  res[is.na(res)] <- 0
  
  if (!.ambig) {
    res <- .fix.ambig(res, .ambig)
  }
  
  if (.norm) {
    if (length(.genes) == 1) {
      res[,-1] <- apply(as.matrix(res[,-1]), 2, function (col) col / sum(col))
    } else {
      res <- res / sum(res)
    }
  }
  
  res
}


#' Perform PCA on segments frequency data.
#' 
#' @aliases pca.segments pca.segments.2D
#' 
#' @description
#' Perform PCA on gene segments frequency data for V- and J-segments and either return pca object or plot the results.
#' 
#' @usage
#' pca.segments(.data, .cast.freq.seg = T, ..., .text = T, .do.plot = T)
#' 
#' pca.segments.2D(.data, .cast.freq.seg = T, ..., .text = T, .do.plot = T)
#' 
#' @param .data Either data.frame or a list of data.frame or a result obtained from the \code{geneUsage} function.
#' @param .cast.freq.seg if T then apply code{geneUsage} to the supplied data.
#' @param ... Further arguments passed to \code{prcomp} or \code{geneUsage}.
#' @param .text If T then plot sample names in the resulting plot.
#' @param .do.plot if T then plot a graphic, else return a pca object.
#' 
#' @return If .do.plot is T than ggplot object; else pca object.
#' 
#' @examples
#' \dontrun{
#' # Load the twins data.
#' data(twb)
#' # Plot a plot of results of PCA on V-segments usage.
#' pca.segments(twb, T, scale. = T)
#' }
pca.segments <- function(.data, .cast.freq.seg = T, ..., .text = T, .do.plot = T){
  if (.cast.freq.seg) { .data <- geneUsage(.data, ...)[,-1] }
  pca.res <- prcomp(t(as.matrix(.data)), ...)
  if (.do.plot) {
    pca.res <- data.frame(PC1 = pca.res$x[,1], PC2 = pca.res$x[,2], Subject = names(.data))
    p <- ggplot() + geom_point(aes(x = PC1, y = PC2, colour = Subject), size = 3, data = pca.res)
    if (.text) {
      p <- p + geom_text(aes(x = PC1, y = PC2, label = Subject), data = pca.res, hjust=.5, vjust=-.3)
    }
    p + theme_linedraw() + guides(size=F) + ggtitle("VJ-usage: Principal Components Analysis") + .colourblind.discrete(length(pca.res$Subject), T)
  } else {
    pca.res
  }
}

pca.segments.2D <- function(.data, .cast.freq.seg = T, ..., .text = T, .do.plot = T){
  if (.cast.freq.seg) { .data <- lapply(geneUsage(.data, ...), function (x) as.vector(x)) }
  pca.res <- prcomp(do.call(rbind, .data), ...)
  if (.do.plot) {
    pca.res <- data.frame(PC1 = pca.res$x[,1], PC2 = pca.res$x[,2], Subject = names(.data))
    p <- ggplot() + geom_point(aes(x = PC1, y = PC2, colour = Subject), size = 3, data = pca.res)
    if (.text) {
      p <- p + geom_text(aes(x = PC1, y = PC2, label = Subject), data = pca.res, hjust=.5, vjust=-.3)
    }
    p <- p + theme_linedraw() + guides(size=F) + ggtitle("VJ-usage: Principal Components Analysis") + .colourblind.discrete(length(pca.res$Subject), T)
  } else {
    pca.res
  }
}