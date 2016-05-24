########## Various plotting functions ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Segment", 'Size', 'Freq', 'Sample', 'V.gene', 'J.gene', '..count..', 'Time.point', 'Proportion', 'Sequence',
                           'Lower', 'Upper', 'Lengths', 'Read.count', 'Var', 'Value', 'Group', 'variable', 'name', 'value', 'Kmers',
                           'Count', 'People', 'First', 'Second', 'Var1', 'Q0.025', 'Q0.975', 'Mean', 'Type', 'Clone.size', 'Q1', 'Q2', 
                           'Symbol', 'Gene', 'Genes', 'Sample', 'label', 'Xrep', 'Yrep'))
}


# red - yellow - green
.ryg.gradient <- function (.min = NA, .max = NA) {
  cs <- c('#66FF00', '#FFFF66', '#FF6633')
  if (!is.na(.min)) {
    scale_fill_gradientn(limits = c(.min, .max), colours = cs, na.value = 'grey60')
  } else {
    scale_fill_gradientn(colours = cs, na.value = 'grey60')
  }
}

# white/orange/yellow - green - blue
# colourblind - friendly
# for fill
.colourblind.gradient <- function (.min = NA, .max = NA, .colour = F) {
  #   cs <- c("#FFFFD9", "#41B6C4", "#225EA8")
#   cs <- c("#FFFFBB", "#41B6C4", "#225EA8")
#   cs <- c("#FFBB00", "#41B6C4", "#225EA8") <- old version
#   cs <- c("#FF4B20", "#FFB433", "#C6EDEC", "#85CFFF", "#0348A6")
  # cs <- c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
  # scale_fill_gradientn(guide='colourbar', colours=c("#0072B2", "#EEEEEE", "#D55E00")
  
  cs <- c(c("#0072B2", "#EEEEEE", "#D55E00"))
  
  if (!is.na(.min)) {
    if (.colour) {
      scale_colour_gradientn(limits = c(.min, .max), guide='colorbar', colours = cs, na.value = 'grey60')
    } else {
      scale_fill_gradientn(limits = c(.min, .max), guide='colorbar', colours = cs, na.value = 'grey60')
    }
  } else {
    if (.colour) {
      scale_colour_gradientn(colours = cs, na.value = 'grey60')
    } else {
      scale_fill_gradientn(colours = cs, na.value = 'grey60')
    }
  }
}

# white/orange/yellow - green - blue
# colourblind - friendly
# for fill and colour
.colourblind.discrete <- function (.n, .colour = F) {
  #   cs <- c("#FFFFD9", "#41B6C4", "#225EA8")
  #   cs <- c("#FFFFBB", "#41B6C4", "#225EA8")
#   cs <- c("#FFBB00", "#41B6C4", "#225EA8") <- old version
  # cs <- c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
  cs <- c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
  if (.colour) {
    scale_colour_manual(values = colorRampPalette(cs)(.n))
  } else {
    scale_fill_manual(values = colorRampPalette(cs)(.n))
  }
}

.colourblind.discrete2 <- function (.n, .colour = F) {
  #   cs <- c("#FFFFD9", "#41B6C4", "#225EA8")
  #   cs <- c("#FFFFBB", "#41B6C4", "#225EA8")
    cs <- c("#FFAB00", "#41B6C4", "#225EA8") # <- old version
  # cs <- c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
  # cs <- c("#FF4B20", "#FFB433", "#0348A6")
  if (.colour) {
    scale_colour_manual(values = colorRampPalette(cs)(.n))
  } else {
    scale_fill_manual(values = colorRampPalette(cs)(.n))
  }
}

# light blues - dark blues
# for fill
.blues.gradient <- function (.min = NA, .max = NA) {
  cs <- c("#F7FBFF", "#9ECAE1", "#2171B5")
  if (!is.na(.min)) {
    scale_fill_gradientn(limits = c(.min, .max), colours = cs, na.value = 'grey60')
  } else {
    scale_fill_gradientn(colours = cs, na.value = 'grey60')
  }
}


#' Plot a histogram of lengths.
#' 
#' @description
#' Plot a histogram of distribution of lengths of CDR3 nucleotide sequences. On y-axis are sum of read counts for each length.
#' 
#' @param .data Data frame with columns 'CDR3.nucleotide.sequence' and 'Read.count' or list with such data frames.
#' @param .ncol If .data is a list, than number of columns in a grid of histograms for each data frame in \code{.data}. Else not used.
#' @param .name Title for this plot.
#' @param .col Name of the column to use in computing the lengths distribution.
#' 
#' @details
#' If \code{.data} is a data frame, than one histogram will be plotted. Is \code{.data} is a list, than grid of histograms
#' will be plotted.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' load('immdata.rda')
#' # Plot one histogram with main title.
#' vis.count.len(immdata[[1]], 'Main title here')
#' # Plot a grid of histograms with 2 columns.
#' vis.count.len(immdata, 2)
#' }
vis.count.len <- function (.data, .ncol = 3, .name = "", .col = 'Read.count') {
  if (has.class(.data, 'list')) {
    return(do.call(grid.arrange, c(lapply(1:length(.data), function (i) vis.count.len(.data[[i]], .col = .col, .name = names(.data)[i])), ncol = .ncol)))
  }
  tmp <- aggregate(as.formula(paste0(.col, " ~ nchar(CDR3.nucleotide.sequence)")), .data, sum)
  names(tmp) <- c('Lengths', "Count")
  ggplot() +
    geom_bar(aes(x = Lengths, y = Count, fill = Count), data = tmp, stat = 'identity', colour = 'black') +
    .colourblind.gradient(min(tmp$Count), max(tmp$Count)) +
    ggtitle(.name) + theme_linedraw()
}


#' Plot a histogram of counts.
#' 
#' @description
#' Plot a histogram of distribution of counts of CDR3 nucleotide sequences. On y-axis are number of counts.
#' 
#' @param .data Cloneset data frame or a list of clonesets.
#' @param .ncol If .data is a list, than number of columns in a grid of histograms for each data frame in \code{.data}. Else not used.
#' @param .name Title for this plot.
#' @param .col Name of the column with counts.
#' 
#' @details
#' If \code{.data} is a data frame, than one histogram will be plotted. Is \code{.data} is a list, than grid of histograms
#' will be plotted.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' load('immdata.rda')
#' # Plot one histogram with main title.
#' vis.number.count(immdata[[1]], 'Main title here')
#' # Plot a grid of histograms with 2 columns.
#' vis.number.count(immdata, 2)
#' }
vis.number.count <- function (.data, .ncol = 3, .name = 'Histogram of clonotypes read counts', .col = "Read.count") {
#   cat('Limits for x-axis set to (0,50). Transform y-axis to sqrt(y).\n')
  
  if (has.class(.data, 'list')) {
    return(do.call(grid.arrange, c(lapply(1:length(.data), function (i) vis.number.count(.data[[i]], .col = .col, .name = names(.data)[i])), ncol = .ncol)))
  }
  
  counts <- data.frame(Count = .data[[.col]])
  
  ggplot() + 
    xlim(min(counts$Count), 300) + 
    ylab('Frequency') +
    geom_histogram(aes(x = Count, fill = ..count..), data = counts, binwidth = 1, colour = 'black') +
    coord_trans(x = 'log10') + scale_y_log10() +
    ggtitle(.name) + 
    .colourblind.gradient() +
    theme_linedraw()
}


#' Heatmap.
#' 
#' @aliases vis.heatmap
#' 
#' @description
#' Plot a heatmap from a matrix or a data.frame
#'
#' @param .data Either a matrix with colnames and rownames specifyed or a data.frame with the first column of
#' strings for row names and other columns stands for values.
#' @param .title Main title of the plot.
#' @param .labs Labs names. Character vector of length 2 (for naming x-axis and y-axis).
#' @param .legend Title for the legend.
#' @param .na.value Replace NAs with this values.
#' @param .text if T then print \code{.data} values at tiles.
#' @param .scientific If T then force show scientific values in the heatmap plot.
#' @param .signif.digits Number of significant digits to show. Default - 4.
#' @param .size.text Size for the text in the cells of the heatmap, 4 by default.
#' @param .no.legend If T than remove the legend from the plot.
#' @param .no.labs If T than remove x / y labels names from the plot.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' # Load your data.
#' load('immdata.rda')
#' # Perform cloneset overlap by amino acid sequences with V-segments.
#' imm.av <- repOverlap(immdata, .seq = 'aa', .vgene = T)
#' # Plot a heatmap.
#' vis.heatmap(imm.av, .title = 'Immdata - (ave)-intersection')
#' }
vis.heatmap <- function (.data, 
                         .title = "Number of shared clonotypes", 
                         .labs = c('Sample', 'Sample'), 
                         .legend = 'Shared clonotypes', 
                         .na.value = NA, 
                         .text = T, 
                         .scientific = FALSE, 
                         .signif.digits = 4,
                         .size.text = 4, 
                         .no.legend = F, 
                         .no.labs = F) {
  if (has.class(.data, 'data.frame')) {
    names <- .data[,1]
    .data <- as.matrix(.data[,-1])
    row.names(.data) <- names
  }

  if (is.null(colnames(.data))) {
    colnames(.data) <- paste0('C', 1:ncol(.data))
  }
  
  if (is.null(row.names(.data))) {
    row.names(.data) <- paste0('C', 1:nrow(.data))
  }
  
  .data[is.na(.data)] <- .na.value
  
  tmp <- as.data.frame(.data)
  tmp$name <- row.names(.data)

  m <- melt(tmp, id.var = c('name'))
  m[,1] <- factor(m[,1], levels = rev(rownames(.data)))
  m[,2] <- factor(m[,2], levels = colnames(.data))
  
  .cg <- .colourblind.gradient(min(m$value), max(m$value))

  m$label <- format(m$value, scientific = .scientific, digits = .signif.digits)
  
  p <- ggplot(m, aes(x = variable, y = name, fill = value))
  p <- p + geom_tile(aes(fill = value), colour = "white")
  if (.text) {
    p <- p + geom_text(aes(fill = value, label = label), size = .size.text)
  }
#   p <- p + geom_text(aes(fill = value, label = value))
#   p <- p + .ryg.gradient(min(m$value), max(m$value))
  p <- p + .cg
  # p <- p + .blues.gradient(min(m$value), max(m$value))
  
  p <- p + ggtitle(.title) + 
    guides(fill = guide_colourbar(title=.legend)) +
    xlab(.labs[1]) + ylab(.labs[2]) + coord_fixed() +
    theme_linedraw() + theme(axis.text.x  = element_text(angle=90)) +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))
  
  if (.no.legend) {
    p <- p + theme(legend.position="none")
  }
  
  if (.no.labs) {
    p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  
  p
}


#' Boxplot for groups of observations.
#' 
#' @aliases vis.group.boxplot
#' 
#' @description
#' Plot boxplots for each group.
#'
#' @param .data Either a matrix with colnames and rownames specifyed or a data frame with the first column of
#' strings for row names and other columns stands for values.
#' @param .groups Named list with character vectors for names of elements for each group. If NA than each
#' member is in the individual group.
#' @param .title Main title of the plot.
#' @param .labs Labs names. Character vector of length 1 (for naming both axis with same name) or 2 (first elements stands for x-axis).
#' @param .rotate.x if T then rotate x-axis.
#' @param .violin If T then plot a violin plot.
#' @param .notch "notch" parameter to the \code{geom_boxplot} ggplo2 function.
#' @param ... Parameters passed to \code{melt}, applied to \code{.data} before plotting in \code{vis.group.boxplot}.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' names(immdata)  # "A1" "A2" "B1" "B2" "C1" "C2"
#' # Plot a boxplot for V-usage for each plot
#' # three boxplots for each group.
#' vis.group.boxplot(freq.Vb(immdata),
#'    list(A = c('A1', 'A2'), B = c('B1', 'B2'), C = c('C1', 'C2')),
#'    c('V segments', 'Frequency')) 
#' 
#' data(twb)
#' ov <- repOverlap(twb)
#' sb <- matrixSubgroups(ov, list(tw1 = c('Subj.A', 'Subj.B'), tw2 = c('Subj.C', 'Subj.D')));
#' vis.group.boxplot(sb)
#' }
vis.group.boxplot <- function (.data, .groups = NA, .labs = c('V genes', 'Frequency'), .title = '', .rotate.x = T, .violin = T, .notch = F, ...) {
  if (has.class(.data, 'data.frame')) {
    .data$Sample <- .data[,1]
    .data <- .data[,c(1,3,2)]
  } else {
    .data <- melt(.data, ...)
  }
  
  colnames(.data) <- c('Var', 'Sample', 'Value')
  .data$Group <- as.character(.data$Sample)
  if (!is.na(.groups)[1]) {
    for (i in 1:length(.groups)) {
      for (name in .groups[[i]]) {
        .data$Group[.data$Sample == name] <- names(.groups)[i]
      }
    }
  }
  
  p <- ggplot() + 
    geom_boxplot(aes(x = Var, y = Value, fill = Group), data = .data, colour = 'black', notch = .notch)
  
  if (.violin) {
    p <- p +geom_violin(aes(x = Var, y = Value, fill = Group), alpha = .2, data = .data)
  }
    
  if (length(.labs) >= 2) {
    p <- p + xlab(.labs[1]) + ylab(.labs[2])
  }
  p <- p + ggtitle(.title) + theme_linedraw()
  if (.rotate.x) {
    p <- p + theme(axis.text.x  = element_text(angle=90))
  }
  p + .colourblind.discrete(length(unique(.data$Group)))
}


#' Histogram of segments usage.
#' 
#' @aliases vis.V.usage vis.J.usage
#' 
#' @description
#' Plot a histogram or a grid of histograms of V- / J-usage.
#' 
#' @param .data Mitcr data frame or a list with mitcr data frames.
#' @param .genes Gene alphabet passed to \link{geneUsage}.
#' @param .main Main title of the plot.
#' @param .ncol Number of columns in a grid of histograms if \code{.data} is a list and \code{.dodge} is F.
#' @param .coord.flip if T then flip coordinates.
#' @param .labs Character vector of length 2 with names for x-axis and y-axis.
#' @param .dodge If \code{.data} is a list, than if this is T plot V-usage for all data frames to the one histogram.
#' @param ... Parameter passed to \code{geneUsage}. By default the function compute V-usage or J-usage for beta chains
#' w/o using read counts and w/ "Other" segments.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' # Load your data.
#' load('immdata.rda')
#' # Compute V-usage statistics.
#' imm1.vs <- geneUsage(immdata[[1]], HUMAN_TRBV)
#' vis.V.usage(immdata, HUMAN_TRBV, .main = 'Immdata V-usage [1]', .dodge = T)
#' # Plot a histogram for one data frame using all gene segment data from V.gene column.
#' vis.V.usage(imm1.vs, NA, .main = 'Immdata V-usage [1]')
#' # Plot a grid of histograms - one histogram for V-usage for each data frame in .data.
#' vis.V.usage(immdata, HUMAN_TRBV, .main = 'Immdata V-usage', .dodge = F, .other = F)
#' }
vis.gene.usage <- function (.data, .genes = NA, .main = "Gene usage", .ncol = 3, .coord.flip = F, .dodge = F, .labs = c("Gene", "Frequency"), ...) {  
  if (!is.na(.genes[1])) {
    res <- geneUsage(.data, .genes, ...)
  } else {
    res <- .data
  }
  
  if (class(res[[2]]) != "factor") {
    res <- melt(res)
    res <- res[1:nrow(res), ] 
    colnames(res) <- c('Gene', 'Sample', 'Freq')
  }
  
  if (length(unique(res$Sample)) > 1) {    
    if (.dodge) {      
      ggplot() + 
        geom_bar(aes(x = Gene, y = Freq, fill = Sample), data = res, stat = 'identity', position = position_dodge(), colour = 'black') +
        theme_linedraw() + 
        theme(axis.text.x = element_text(angle=90)) + 
        .colourblind.discrete(length(unique(res$Sample))) +
        scale_y_continuous(expand = c(.02,0))
    } else {
      res <- split(res, res$Sample)
      ps <- lapply(1:length(res), function (i) {
        vis.gene.usage(res[[i]], NA, names(res)[i], 0, .coord.flip, .labs = .labs, ...) 
      })      
      do.call(grid.arrange, c(ps, ncol = .ncol, top = .main) )
    }
    
  }  
  else {        
    p <- ggplot() + 
      geom_bar(aes(x = Gene, y = Freq, fill = Freq), data = res, stat = 'identity', colour = 'black')
    
    if (.coord.flip) { p <- p + coord_flip() }
    
    p + theme_linedraw() + 
      theme(axis.text.x = element_text(angle=90)) + 
      ggtitle(.main) + 
      .colourblind.gradient() +
      scale_y_continuous(expand = c(.02,0)) + 
      xlab(.labs[1]) + ylab(.labs[2])
  }
}


#' PCA result visualisation
#' 
#' @description
#' Plot the given pca results with colour divided by the given groups.
#' 
#' @param .data Result from prcomp() function or a data frame with two columns 'First' and 'Second'
#' stands for the first PC and the second PC.
#' @param .groups List with names for groups and indices of the group members. If NA than each
#' member is in the individual group.
#' 
#' @return ggplot object.
vis.pca <- function (.data, .groups = NA) {
  if (has.class(.data, 'data.frame')) {
    dnames <- row.names(.data)
    .data <- data.frame(First = .data[,1], Second = .data[,2], Sample = row.names(.data),
                        Group = rep('group0', times = length(.data[,2])), stringsAsFactors=F)
  } else {
    dnames <- row.names(.data$x)
    .data <- data.frame(First = .data$x[,1], Second = .data$x[,2], Sample = row.names(.data$x),
                        Group = rep('group0', times = length(.data$x[,2])), stringsAsFactors=F)
  }
  
  if (is.na(.groups[1])) { 
    .groups <- lapply(1:nrow(.data), function (i) i)
    names(.groups) <- dnames
  }
  for (i in 1:length(.groups)) {
    for (j in 1:length(.groups[[i]])) {
      .data$Group[.groups[[i]][j]] <- names(.groups)[i]
    }
  }
  
  ggplot() + 
    geom_point(aes(x = First, y = Second, colour = Group), size = 3, data = .data) + 
    geom_text(aes(x = First, y = Second, label = Sample, colour = Group), data = .data, hjust=0, vjust=0) +
    theme_linedraw() +
    .colourblind.discrete2(length(.groups), T)
}


#' Radar-like / spider-like plots.
#' 
#' @description
#' Plot a grid of radar(-like) plots for visualising a distance among objects.
#' 
#' @param .data Square data frame or matrix with row names and col names stands for objects and values for distances.
#' @param .ncol Number of columns in the grid.
#' @param .expand Interger vector of length 2, for \code{scale_y_continous(expand = .expand)} function.
#' 
#' @seealso \link{repOverlap}, \link{js.div}
#' 
#' @examples
#' \dontrun{
#' load('immdata.rda')
#' # Compute Jensen-Shannon divergence among V-usage of repertoires.
#' imm.js <- js.div.seg(immdata, .verbose = F)
#' # Plot it.
#' vis.radarlike(imm.js)
#' }
vis.radarlike <- function (.data, .ncol = 3, .expand = c(.25, 0)) {
  step = ncol(.data)
  data.names <- colnames(.data)
  .data <- as.data.frame(melt(.data))
  .data[is.na(.data[,3]),3] <- 0
  ps <- lapply(seq(1, nrow(.data), step), function (l) {
    ggplot(.data[l:(l+step-1),], aes(x = Var1, y = value, fill = Var1)) + 
      geom_bar(colour = 'black', stat = 'identity') + 
      coord_polar() + 
      ggtitle(names(.data)[l]) +
      scale_y_continuous(expand = .expand) + 
      guides(fill = guide_legend(title="Sample")) +
      theme_linedraw() + xlab('') + ylab('')
    })
  for (i in 1:length(data.names)) {
    ps[[i]] <- ps[[i]] + ggtitle(data.names[i]) + .colourblind.discrete(length(data.names))
  }
  
  do.call(grid.arrange, c(ps, ncol = .ncol))
}


#' Visualisation of top clones proportions.
#' 
#' @description
#' Visualisation of proportion of the top clones.
#' 
#' @param .data Data frame with clones.
#' @param .head Integer vector of clones for the \code{.head} parameter for the \code{top.proportion} function.
#' @param .col Parameter \code{.col} for the \code{top.proportion} function.
#' 
#' @seealso \code{top.proportion}
#' 
#' @examples
#' \dontrun{
#' vis.top.proportions(immdata)
#' }
vis.top.proportions <- function (.data, .head = c(10, 100, 1000, 10000, 30000, 100000, 300000, 1000000), .col = "Read.count") {
  if (has.class(.data, 'data.frame')) {
    .data <- list(Sample = .data)
  }
  
  res <- sapply(.head, function (h) top.proportion(.data, h, .col))
  tmp <- res
  if (is.null(dim(tmp))) {
    tmp <- t(as.matrix(tmp))
    res <- t(as.matrix(res))
  }
  for (i in 2:ncol(res)) {
    tmp[,i] <- res[,i] - res[,i-1]
  }
  res <- tmp
  colnames(res) <- paste0('[', c(1, .head[-length(.head)] + 1), ':', .head, ')')
  res <- as.data.frame(res)
  res$People <- factor(row.names(res), levels = row.names(res))
  res <- melt(res)
  #   res$variable <- factor(as.character(res$variable), labels = paste0('[', c(1, .head[-length(.head)] + 1), ':', .head, ')'), ordered = T)
  ggplot() + geom_bar(aes(x = People, y = value, fill = variable), data = res, stat = 'identity', position = 'stack', colour = 'black')+ 
    theme_linedraw()  + 
    theme(axis.text.x  = element_text(angle=90)) +
    ylab("Clonal proportion") + 
    xlab("Sample") + 
    ggtitle("Summary proportion of the top N clones")  + 
    guides(fill = guide_legend("Top N clones")) + .colourblind.discrete(length(.head))
#     scale_y_continuous(expand = c(0, 0))
}


#' Rarefaction statistics visualisation.
#' 
#' @description
#' Plot a line with mean unique clones.
#' 
#' @param .muc.res Output from the \code{muc} function.
#' @param .groups List with names for groups and names of the group members. If NULL than each
#' member is in the individual group.
#' @param .log if T then log-scale the y axis.
#' 
#' @seealso \link{rarefaction}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' names(twb)  # "Subj.A" "Subj.B" "Subj.C" "Subj.D"
#' twb.rar <- rarefaction(twb, .col = "Read.count")
#' vis.rarefaction(twb.rar, list(A = c("Subj.A", "Subj.B"), B = c("Subj.C", "Subj.D")))
#' }
vis.rarefaction <- function (.muc.res, .groups = NULL, .log = F) {
  .muc.res$Group <- .muc.res$People
  
  if (!is.null(.groups)) { 
    for (i in 1:length(.groups)) {
      for (j in 1:length(.groups[[i]])) {
        .muc.res$Group[.muc.res$People == .groups[[i]][j] ] <- names(.groups)[i]
      }
    }
  }
  
  .muc.res$Type <- factor(.muc.res$Type, levels = c('interpolation', 'extrapolation'), ordered = T)
  
  p <- ggplot() + 
#     geom_point(aes(x = Size, y = Mean, colour = Group), data = .muc.res, size = 2) + 
    geom_line(aes(x = Size, y = Mean, colour = Group, Group = People, linetype = Type), data = .muc.res) + 
#     geom_errorbar(aes(x = Size, y = Mean, ymin = Q0.025, ymax = Q0.975, colour = Group), data = .muc.res) +
    xlab('Sample size') + ylab('Clones') + ggtitle("Rarefaction analysis") +
    theme_linedraw() + .colourblind.discrete(length(unique(.muc.res$Group)), T)
  
  for (subj in unique(.muc.res$People)) {
    tmp <- tail(.muc.res[.muc.res$People == subj, ], 1)
    p <- p + geom_text(aes(x = Size, y = Mean, label = People), data = tmp, hjust=1, vjust=1)
  }
  
  if (.log) {
    p <- p + scale_x_log10()
  }
  p
}


#' Plot of the most frequent kmers.
#' 
#' @description
#' Plot a distribution (bar plot) of the most frequent kmers in a data.
#' 
#' @param .kmers Data frame with two columns "Kmers" and "Count" or a list with such data frames. See Examples.
#' @param .head Number of the most frequent kmers to choose for plotting from each data frame.
#' @param .position Character vector of length 1. Position of bars for each kmers. Value for the \code{ggplot2} argument \code{position}.
#' 
#' @seealso \code{get.kmers}
#' 
#' @examples
#' \dontrun{
#' # Load necessary data and package.
#' library(gridExtra)
#' load('immdata.rda')
#' # Get 5-mers.
#' imm.km <- get.kmers(immdata)
#' # Plots for kmer proportions in each data frame in immdata.
#' p1 <- vis.kmer.histogran(imm.km, .position = 'stack')
#' p2 <- vis.kmer.histogran(imm.km, .position = 'fill')
#' grid.arrange(p1, p2)
#' }
vis.kmer.histogram <- function (.kmers, 
                                .head = 100, 
                                .position = c('stack', 'dodge', 'fill')) {
  kmers.df <- data.frame(Kmers = '')
  for (i in 2:ncol(.kmers)) {
    kmers.df <- merge(head(.kmers[order(.kmers[, i], decreasing = T), c(1,i)], .head), kmers.df, all = T)
  }
  kmers.df[is.na(kmers.df)] <- 0
  kmers.df <- melt(kmers.df[-1,])
  names(kmers.df) <- c('Kmers', 'People', 'Count')
  p <- ggplot() + geom_bar(aes(x = Kmers, y = Count, fill = People), data = kmers.df, stat = 'identity', position = .position[1]) + theme_linedraw()
  if (.position[1] == 'stack' || .position[1] == 'dodge') {
    p <- p + ylab('Count') + theme(axis.text.x  = element_text(angle=90))
  } else {
    p <- p + ylab('Proportions') + theme(axis.text.x  = element_text(angle=90))
  }
  p + scale_y_continuous(expand = c(0, 0)) + .colourblind.discrete2(length(unique(kmers.df$People)))
}


#' Visualise clonal dynamics among time points.
#' 
#' @description
#' Visualise clonal dynamics (i.e., changes in frequency or count) with error bars of given
#' clones among time points.
#' 
#' @param .changed Result from the \code{find.clonotypes} function, i.e. data frame with first
#' columns with sequences (nucleotide or amino acid) and other columns are columns with frequency / count
#' for each time point for each clone.
#' @param .lower Similar to .changed but values are lower bound for clonal count / frequency.
#' @param .upper Similar to .changed but values are upper bound for clonal count / frequency.
#' @param .log if T then log-scale y-axis.
#' 
#' @return ggplot object.
vis.clonal.dynamics <- function (.changed, .lower, .upper, .log = T) {
  .changed <- melt(.changed, id.vars = names(.changed)[1])
  .lower <- melt(.lower, id.vars = names(.changed)[1])
  .upper <- melt(.upper, id.vars = names(.changed)[1])
  names(.changed) <- c('Sequence', 'Time.point', 'Proportion')
  d <- cbind(.changed, Lower = .lower[,3], Upper = .upper[,3])
  p <- ggplot() + geom_line(aes(x = Time.point, y = Proportion, colour = Sequence, group = Sequence), data = d) +
    geom_errorbar(aes(x = Time.point, y = Proportion, colour = Sequence, ymin = Lower, ymax = Upper), data = d, width = .25) +
    theme_linedraw() + theme(axis.text.x  = element_text(angle=90)) +
    .colourblind.discrete(length(unique(.changed$Sequence)), .colour = T)
  if (.log) {
    p <- p + scale_y_log10()
  }
  p
}


#' Visualise occupied by clones homeostatic space among Samples or groups.
#' 
#' @description
#' Visualise which clones how much space occupy.
#' 
#' @param .clonal.space.data Data from the \code{fclonal.space.homeostasis} function.
#' @param .groups List of named character vector with names of Samples 
#' in \code{.clonal.space.data} for grouping them together.
#' 
#' @seealso \link{clonal.space.homeostasis}
#' 
#' @return ggplot object.
vis.clonal.space <- function (.clonal.space.data, .groups = NULL) {
  melted <- melt(.clonal.space.data)
  colnames(melted) <- c('Sample', 'Clone.size', 'Proportion')
  melted$Sample <- as.character(melted$Sample)
  melted$Proportion <- as.numeric(as.character(melted$Proportion))
  melted$Group <- melted$Sample
  
  if (!is.null(.groups)) { 
    for (i in 1:length(.groups)) {
      for (j in 1:length(.groups[[i]])) {
        melted$Group[melted$Sample == .groups[[i]][j] ] <- names(.groups)[i]
      }
    }
    
    perc <- melt(tapply(melted$Proportion, list(melted$Group, melted$Clone.size), function (x) c(quantile(x, probs = .25), mean(x), quantile(x, probs = .75))))
    return(perc)
    perc <- data.frame(row.names(perc), perc, stringsAsFactors = F)
    colnames(perc) <- c('Group', 'Q1', 'Mean', 'Q2')
    
    p <- ggplot() +
      geom_bar(aes(x = Group, y = Mean, fill = Clone.size), data = melted, colour = 'black', stat = 'identity') +
      geom_errorbar(aes(x = Group, ymin = Q1, ymax = Q2), data = melted, colour = 'black') +
      xlab("Sample")
  } else {
    p <- ggplot() +
      geom_bar(aes(x = Group, y = Proportion, fill = Clone.size), data = melted, colour = 'black', stat = 'identity', position = 'stack') +
      xlab("Sample")
      
  }
    
  p + theme_linedraw() + 
    theme(axis.text.x = element_text(angle=90)) + ylab("Occupied homeostatic space, proportion") + 
    ggtitle("Clonal space homeostasis") + 
    guides(fill = guide_legend("Clone size")) + .colourblind.discrete(length(unique(melted$Clone.size))) +
    scale_y_continuous(expand = c(.01, .01)) + scale_x_discrete(expand = c(.02, .02))
}


#' Logo - plots for amino acid and nucletide profiles.
#' 
#' @description
#' Plot logo-like graphs for visualising of nucleotide or amino acid motif sequences / profiles.
#' 
#' @param .data Output from the \code{kmer.profile} function.
#' @param .replace.zero.with.na if T then replace all zeros with NAs, therefore letters with
#' zero frequency wont appear at the plot.
#' @param .jitter.width,.jitter.height,.dodge.width Parameters to \code{position_jitterdodge}
#' for aligning text labels of letters.
#' 
#' @return ggplot2 object
#' 
#' @examples
#' \dontrun{
#' d <- kmer.profile(c('CASLL', 'CASSQ', 'CASGL'))
#' vis.logo(d)
#' }
vis.logo <- function (.data, .replace.zero.with.na = T, .jitter.width = .01, .jitter.height = .01, .dodge.width = .15) {
  .data <- melt(.data)
  if (.replace.zero.with.na) {
    .data$value[.data$value == 0] <- NA
  }
  ggplot(aes(x = variable, y = value, fill = Symbol, colour = Symbol), data = .data) + 
    geom_point(colour = 'black') + 
    geom_text(aes(label = Symbol), size = 5, 
              position = position_jitterdodge(jitter.width = .jitter.width, 
                                              jitter.height = .jitter.height, 
                                              dodge.width = .dodge.width)) +
    xlab("Position") + ylab("Proportion") +
    theme_linedraw()
}


#' Visualisation of shared clonotypes occurrences among repertoires.
#' 
#' @description 
#' Visualise counts or proportions of shared clonotypes among repertoires.
#' 
#' @param .shared.rep Shared repertoires, as from \link{shared.repertoire} function.
#' @param .x.rep Which repertoire show on x-axis. Either a name or an index of a repertoire 
#' in the \code{.shared.rep} or NA to choose all repertoires.
#' @param .y.rep Which repertoire show on y-axis. Either a name or an index of a repertoire 
#' in the \code{.shared.rep} or NA to choose all repertoires.
#' @param .title Main title of the plot.
#' @param .ncol Number of columns in the resulting plot.
#' @param .point.size.modif Modify this to correct sizes of points.
#' 
#' @return ggplot2 object or plot
#' 
#' @seealso \link{shared.repertoire}
#' 
#' @examples 
#' \dontrun{
#' data(twb)
#' # Show shared nucleotide clonotypes of all possible pairs 
#' # using the Read.proportion column
#' twb.sh <- shared.repertoire(twb, "n0rp")
#' vis.shared.clonotypes(twb.sh, .ncol = 4)
#' 
#' # Show shared amino acid + Vseg clonotypes of pairs 
#' # including the Subj.A (the first one) using
#' # the Read.count column.
#' twb.sh <- shared.repertoire(twb, "avrc")
#' vis.shared.clonotypes(twb.sh, 1, NA, .ncol = 4)
#' # same, just another order of axis
#' vis.shared.clonotypes(twb.sh, NA, 1, .ncol = 4)
#' 
#' # Show shared nucleotide clonotypes of Subj.A (the first one)
#' # Subj.B (the second one) using the Read.proportion column.
#' twb.sh <- shared.repertoire(twb, "n0rp")
#' vis.shared.clonotypes(twb.sh, 1, 2)
#' 
#' # Show the same plot, but with much larget points.
#' vis.shared.clonotypes(twb.sh, 1, 2, .point.size.modif = 3)
#' }
vis.shared.clonotypes <- function (.shared.rep, .x.rep = NA, .y.rep = NA, 
                                   .title = "Shared clonotypes", .ncol = 3, 
                                   .point.size.modif = 1) {
  mat <- shared.matrix(.shared.rep)
  
  if (is.na(.x.rep) && is.na(.y.rep)) {
    ps <- list()
    for (i in 1:ncol(mat)) {
      for (j in 1:ncol(mat)) {
        ps <- c(ps, list(vis.shared.clonotypes(.shared.rep, i, j, '')))
      }
    }
    do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
  } else if (is.na(.x.rep)) {
    ps <- lapply(1:ncol(mat), function (i) { 
      vis.shared.clonotypes(.shared.rep, i, .y.rep, '') 
      })
    do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
  } else if (is.na(.y.rep)) {
    ps <- lapply(1:ncol(mat), function (j) { 
      vis.shared.clonotypes(.shared.rep, .x.rep, j, '') 
    })
    do.call(grid.arrange, c(ps, ncol = .ncol, top = .title))
  } else {
    if (!is.character(.x.rep)) { .x.rep <- colnames(mat)[.x.rep] }
    if (!is.character(.y.rep)) { .y.rep <- colnames(mat)[.y.rep] }
    
    df <- data.frame(cbind(mat[, .x.rep], mat[, .y.rep]))
    df <- df[!is.na(df[,1]) & !is.na(df[,2]), ]
    freq <- log10(sqrt(as.numeric(df[, 1]) * df[, 2])) / 2
    names(df) <- c("Xrep", "Yrep")
    
    pnt.cols <- log(df[, 1] / df[, 2])
    suppressWarnings(pnt.cols[pnt.cols > 0] <- pnt.cols[pnt.cols > 0] / max(pnt.cols[pnt.cols > 0]))
    suppressWarnings(pnt.cols[pnt.cols < 0] <- -pnt.cols[pnt.cols < 0] / min(pnt.cols[pnt.cols < 0]))
    
    mat.lims <- c(min(as.matrix(df)), max(as.matrix(df)))
    
    ggplot() + 
      geom_point(aes(x = Xrep, y = Yrep, size = freq, fill = pnt.cols), data = df, shape=21) + 
      scale_radius(range = c(.point.size.modif, .point.size.modif * 6)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      theme_linedraw() + 
      .colourblind.gradient(min(pnt.cols), max(pnt.cols)) +
      scale_x_log10() + scale_y_log10() + theme(legend.position="none") +
      coord_fixed(xlim = mat.lims, ylim = mat.lims) +
      xlab(.x.rep) + ylab(.y.rep) + ggtitle(.title)
  }
}


# vis.hill.numbers <- function (.hill.nums, .groups = NA) {
#   
# }