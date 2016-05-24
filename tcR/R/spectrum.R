########## Spectratyping ##########


#' Spectratype plot.
#' 
#' @description
#' General function for making a spectratyping plots.
#' 
#' @param .data mitcr List with data frames.
#' @param .column Character vector with name of the column with numeric characteristic.
#' @param .by.alphabet Either 'Va/b', 'Ja/b' or an alphabet.
#' @param .do.legend if T then plot a legend.
#' @param .draw.only.legend if T then plot only a legend without plots.
#' @param .legend.ncol Number of columns in the legend. If -1 than function will try to predict number of columns.
#' @param .nrow Rows of grid of plots.
#' @param .by.col Character vector with name of the column by which divide the given data frames.
#' @param .sum.col Which column use for sum.
#' @param .other if T then include in the result plot values which isn't in the given alphabet.
#' @param .log if T then scale y-axis by log10.
#' @param .verbose if T then print messages about state of the process.
#' 
#' @details
#' For each element in \code{.data} do: for each factor in \code{.by.col} which is in \code{.by.alphabet}, compute histogram of \code{.column}
#' and then plot stacked histogram of distributions of \code{.by.col} for each factor.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' # Spectratyping of distribution of length of CDR3 nucleotide sequences
#' # by V-beta-segments.
#' immdata <- lapply(immdata, function (x) { 
#'              x$Length <- nchar(x$CDR3.nucleotide.sequence)
#'              x
#'            } )
#' spectratyping(immdata, 'Length', 'Vb')
#' # Spectratyping of distribution of Total insertions
#' # by J-beta-segments.
#' spectratyping(immdata, 'Total.insertions', 'Jb')
#' }
spectratyping <- function (.data, .column = 'VD.insertions', .by.alphabet = 'Vb', 
                           .do.legend = T, .draw.only.legend = F,  .legend.ncol = -1,
                           .nrow = 2, .by.col = '', .sum.col = 'Percentage',
                           .other = F, .log = F, .verbose = T) {  
  .by <- .by.alphabet
  if (.by.col == '') {
    if (.by == 'Vb') {
      .by <- HUMAN_TRBV_MITCR
      .by.col <- 'V.gene'
    } else if (.by == 'Va') {
      .by <- HUMAN_TRAV
      .by.col <- 'V.gene'
    } else if (.by == 'Ja') {
      .by <- HUMAN_TRAJ
      .by.col <- 'J.gene'
    } else {
      .by <- HUMAN_TRBJ
      .by.col <- 'J.gene'
    }
  }
  
  if (.legend.ncol == -1) {
    .legend.ncol <- length(.by) %/% 23
    if (.legend.ncol < 1) { .legend.ncol <- 1 }
  }
  
  twb.vs.ins <- list()
  tmp.data <- .data
  
  for (i in 1:length(tmp.data)) {
    if (.verbose) cat(names(tmp.data)[i], '\n')
    
    subdata <- tmp.data[[i]]
    all.len <- unique(subdata[[.column]])
    all.len <- all.len[all.len > -1]
    if (.other) {
      res.table <- matrix(0, length(all.len), length(.by) + 1)
    } else {
      res.table <- matrix(0, length(all.len), length(.by))
    }
    if (.verbose) pb <- set.pb(length(.by) * length(all.len))
    for (k in 1:length(.by)) {
      tmp <- subdata[subdata[[.by.col]] == .by[k],]
      for (j in 1:length(all.len)) {
        res.table[j,k] <- sum(tmp[[.sum.col]][tmp[[.column]] == all.len[j]])
        if (.verbose) add.pb(pb)
      }
    }
    if (.other) {
      tmp <- subdata[!(subdata[[.by.col]] %in% .by),]
      for (j in 1:length(all.len)) {
        res.table[j, length(.by) + 1] <- sum(tmp[[.sum.col]][tmp[[.column]] == all.len[j]])
        if (.verbose) add.pb(pb)
      }
    }
    if (.verbose) close(pb)
    res.table[res.table == 0] <- NA
    res.table <- as.data.frame(res.table)
    res.table <- cbind(all.len, res.table)
    if (.other) {
      names(res.table) <- c('Len', .by, 'Other')
    } else {
      names(res.table) <- c('Len', .by)
    }
    twb.vs.ins[[i]] <- melt(res.table, id.vars ='Len')
    names(twb.vs.ins[[i]]) <- c(.column, .by.col, .sum.col)
  }
  colorsr <- rainbow(length(.by) + 1, s=.6, v=.9)[sample(1:length(.by),length(.by))]
  ps <- lapply(1:length(tmp.data), function (x) {
    p <- ggplot(data = twb.vs.ins[[x]], aes_string(x = .column, y = .sum.col, fill = .by.col)) +
      ggtitle(names(tmp.data)[x]) +
      geom_histogram(stat = 'identity') +
      theme(legend.position = 'none') +
      scale_fill_manual(values=colorsr)
    if (.log) {p <- p +  scale_y_log10()}
    p
    })
  if (.do.legend) {
    leg <- gtable_filter(ggplot_gtable(ggplot_build(ggplot(data = twb.vs.ins[[1]], aes_string(x = .column, y = .sum.col, fill = .by.col)) +
                                                      ggtitle(names(tmp.data)[1]) + geom_histogram(stat = 'identity') +
                                                      guides(fill=guide_legend(ncol=.legend.ncol)) + scale_fill_manual(values=colorsr))),
                         "guide-box")
    grid.arrange(do.call(arrangeGrob, c(ps, nrow = .nrow)), leg, widths=unit.c(unit(1, "npc") - leg$width, leg$width), nrow = 1, top = paste0(.column, ', ', .by.col))
  } else if (!.draw.only.legend) {
    grid.arrange(do.call(arrangeGrob, c(ps, nrow = .nrow)), top = paste0(.column, ', ', .by.col))
  } else {
    grid.arrange(leg)
  }
}

# tmp <- lapply(bek[c(1:6,10:15)], function(x) { x$Length <- nchar(x$CDR3.amino.acid.sequence); x } )
# spectratyping(tmp, 'Length', 'Vb')
# spectratyping(tmp, 'Length', 'Jb')
# spectratyping(bek.del[c(1:6, 10:15)], 'VD.insertions', 'Vb', .nrow = 4)
# spectratyping(bek.del[c(1:6, 10:15)], 'DJ.insertions', 'Jb', .nrow = 4)
# spectratyping(bek.del[c(1:6, 10:15)], 'VD.deletions', 'Vb', .nrow = 2)
# spectratyping(bek.del[c(1:6, 10:15)], 'DJ.deletions', 'Jb', .nrow=2)
# spectratyping(bek.del[1], 'DJ.deletions', .nrow=2, .other = T)

# EXAMPLE HERE
# spectratyping(bek.del[1], 'DJ.deletions', .by.alphabet = head(bek.del[[1]]$CDR3.amino.acid.sequence, 50), .by.col = 'CDR3.amino.acid.sequence',  .nrow=2)
# EXAMPLE WITH READ.COUNT
# spectratyping(bek.del[1:6], 'DJ.deletions', .by.alphabet = head(bek.del[[1]]$CDR3.amino.acid.sequence, 50), .by.col = 'CDR3.amino.acid.sequence',  .nrow=2, .sum.col='Read.count')